"""Exploratory hydration geometry analysis for PA6 deacylation trajectories."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import netcdf_file

from nyl12_covalent.deacylation_water_geometry import (
    ASP_RESIDUES,
    ATTACK_C,
    BASE,
    CARBONYL_O,
    LEAVING_N,
    NALPHA,
    OG1,
    _extension_sources,
    _read_frames,
    _topology,
)
from nyl12_covalent.deacylation_water_reporting import write_reports
from nyl12_covalent.deacylation_water_measurements import measure_frames


def analyze(output: Path, base: Path = BASE) -> dict[str, Any]:
    """Analyze saved endpoint, heating, and completed equilibration trajectories."""
    top_path = base / "inputs/system.prmtop"
    endpoint = base / "inputs/endpoint.ncrst"
    trajectory_paths = [
        (base / "heat/trajectory.nc", "heat"),
        (base / "equil/trajectory.nc", "equil_300K"),
    ] + _extension_sources(base)
    topology, waters = _topology(top_path)
    # Validate mapped residues/names before measuring.
    for idx, name in [
        (ATTACK_C, "C"),
        (CARBONYL_O, "O"),
        (NALPHA, "N"),
        (OG1, "OG1"),
        (LEAVING_N, "N"),
    ]:
        if topology.atoms[idx - 1].name != name:
            raise ValueError(f"Atom mapping mismatch at {idx}: expected {name}")
    for residue_idx, oxygens in ASP_RESIDUES.items():
        residue = topology.residues[residue_idx - 1]
        if residue.name != "ASP" or any(
            topology.atoms[i - 1].name not in {"OD1", "OD2"} for i in oxygens
        ):
            raise ValueError(f"Asp mapping mismatch at residue {residue_idx}")
    nalpha_h = sorted(
        atom.idx
        for atom in topology.atoms[NALPHA - 1].bond_partners
        if atom.element == 1 or atom.name.startswith("H")
    )
    nleave_h = sorted(
        atom.idx
        for atom in topology.atoms[LEAVING_N - 1].bond_partners
        if atom.element == 1 or atom.name.startswith("H")
    )
    with netcdf_file(endpoint, "r", mmap=False) as nc:
        t0 = float(np.asarray(nc.variables["time"].data))
        xyz0 = np.asarray(nc.variables["coordinates"][:])
        box0 = np.concatenate(
            (
                np.asarray(nc.variables["cell_lengths"][:]),
                np.asarray(nc.variables["cell_angles"][:]),
            )
        )
    if (
        not math.isfinite(t0)
        or not np.isfinite(xyz0).all()
        or not np.isfinite(box0).all()
    ):
        raise ValueError(
            "Endpoint restart contains nonfinite time, coordinates, or box"
        )
    frames = [{"stage": "endpoint", "time_ps": t0, "xyz": xyz0, "box": box0}]
    for path, stage in trajectory_paths:
        frames.extend(_read_frames(path, stage))
    all_times = np.asarray([frame["time_ps"] for frame in frames])
    if not np.isfinite(all_times).all() or (
        len(all_times) > 1 and not np.all(np.diff(all_times) > 0)
    ):
        raise ValueError(
            "Combined endpoint and accepted trajectory times must be strictly increasing"
        )
    if any(frame["xyz"].shape != (len(topology.atoms), 3) for frame in frames):
        raise ValueError("Trajectory coordinate atom count does not match topology")
    output.mkdir(parents=True, exist_ok=True)
    nearest_rows, minima = measure_frames(
        frames, waters, topology, nalpha_h, nleave_h, output
    )
    frame_counts: dict[str, int] = {}
    for frame in frames:
        frame_counts[frame["stage"]] = frame_counts.get(frame["stage"], 0) + 1
    for item in minima.values():
        item.pop("_min_dist", None)
        item["near_attack_fraction_by_stage"] = {
            stage: count / frame_counts[stage]
            for stage, count in item["stages"].items()
        }
    ranked = sorted(
        waters,
        key=lambda w: (
            minima[w["o"]]["min_attack_distance_A"],
            -minima[w["o"]]["near_attack_frame_count"],
        ),
    )
    water_rows = []
    for w in ranked:
        m = minima[w["o"]]
        water_rows.append(
            {
                "water_residue_index_1based": w["residue"],
                "water_residue_name": w["name"],
                "oxygen_atom_index_1based": w["o"] + 1,
                "hydrogen_atom_indices_1based": [x + 1 for x in w["h"]],
                **m,
            }
        )
    summary = {
        "schema_version": 1,
        "source_base": str(base),
        "topology": str(top_path),
        "endpoint": str(endpoint),
        "trajectory_sources": [
            {"path": str(p), "stage": s} for p, s in trajectory_paths
        ],
        "frame_count_including_endpoint": len(frames),
        "frame_count_by_stage": frame_counts,
        "pending_extension_segments": [
            {
                "segment": path.parent.name,
                "status": json.loads(path.read_text()).get("status"),
            }
            for path in sorted(base.glob("extend_*/status.json"))
            if json.loads(path.read_text()).get("status") != "accepted"
        ],
        "time_units": "ps",
        "endpoint_time_ps": t0,
        "criteria": {
            "near_attack_O_C_A_max": 3.5,
            "attack_angle_deg": [95, 125],
            "hbond_DA_A_max": 3.5,
            "hbond_DHA_deg_min": 135,
            "interpretation": "exploratory geometry annotation only; does not establish reactivity or proton transfer; water-H...Nalpha is descriptive because Nalpha is protonated",
        },
        "validated_mapping": {
            "attack_C": ATTACK_C,
            "carbonyl_O": CARBONYL_O,
            "Nalpha": NALPHA,
            "OG1": OG1,
            "leaving_N": LEAVING_N,
            "Asp_residue_indices": list(ASP_RESIDUES),
            "Nalpha_H_atom_indices_1based": [index + 1 for index in nalpha_h],
            "Nleaving_H_atom_indices_1based": [index + 1 for index in nleave_h],
        },
        "top5_waters": water_rows[:5],
        "water_count": len(waters),
    }
    accepted_dirs = [
        base / stage for _, stage in trajectory_paths if stage.startswith("extend_")
    ]
    source_paths = [top_path, endpoint, *[path for path, _ in trajectory_paths]]
    source_paths += list(base.glob("extend_*/status.json"))
    source_paths += [
        directory / name
        for directory in accepted_dirs
        for name in ("status.json", "scientific_acceptance.json")
    ]
    code_paths = [
        Path(__file__),
        Path(__file__).with_name("deacylation_water_geometry.py"),
        Path(__file__).with_name("deacylation_water_measurements.py"),
        Path(__file__).with_name("deacylation_water_reporting.py"),
        Path(__file__).parents[2] / "scripts/analyze_deacylation_hydration.py",
    ]
    write_reports(
        output,
        summary,
        water_rows,
        nearest_rows,
        minima,
        len(frames),
        code_paths,
        source_paths,
    )
    return summary
