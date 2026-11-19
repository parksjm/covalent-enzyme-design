"""Screen joint-eligible WT MD frames for native Thr233 alignment."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import shlex
import subprocess
from typing import Any

import numpy as np
import parmed
from parmed.amber import NetCDFTraj

from nyl12_covalent.catalytic_alignment import (
    BACKBONE,
    FIT_RESIDUES,
    SITE_RESIDUES,
    _cif_atoms,
    _topology_atoms,
)
from nyl12_covalent.catalytic_alignment_geometry import fit_coordinates, state_metrics

LOGGER = logging.getLogger(__name__)
NATIVE_THR_LIMITS = {
    "og1_displacement_angstrom": 0.75,
    "n_displacement_angstrom": 1.25,
    "absolute_psi_delta_degrees": 35.0,
    "absolute_chi1_delta_degrees": 20.0,
}


def _source_record(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path.resolve()), "sha256": digest.hexdigest()}


def eligible_frames(path: Path, replica: str) -> list[dict[str, Any]]:
    """Read all joint Tier A/B and native-torsion frames from a prior screen."""
    required = {
        "frame",
        "time_ns",
        "tier",
        "native_torsions",
        "psi_delta",
        "chi1_delta",
    }
    selected: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Screen CSV is missing columns: {sorted(required)}")
        for ordinal, row in enumerate(reader, start=1):
            frame_value = float(row["frame"])
            if not np.isfinite(frame_value) or not frame_value.is_integer():
                raise ValueError(
                    f"Invalid frame index at row {ordinal}: {row['frame']}"
                )
            frame = int(frame_value)
            if frame != ordinal:
                raise ValueError(f"Screen frame indices are not sequential in {path}")
            time_ns = float(row["time_ns"])
            if not np.isfinite(time_ns) or time_ns < 0:
                raise ValueError(f"Invalid time at frame {frame}")
            if row["tier"] not in {"A", "B"}:
                continue
            if row["native_torsions"].lower() not in {"true", "false"}:
                raise ValueError(f"Invalid native_torsions value at frame {frame}")
            psi_delta, chi1_delta = float(row["psi_delta"]), float(row["chi1_delta"])
            native = abs(psi_delta) <= 35.0 and abs(chi1_delta) <= 20.0
            if not np.isfinite([psi_delta, chi1_delta]).all():
                raise ValueError(f"Nonfinite torsion delta at frame {frame}")
            if native != (row["native_torsions"].lower() == "true"):
                raise ValueError(
                    f"Native-torsion flag disagrees with deltas at frame {frame}"
                )
            if native:
                selected.append(
                    {
                        "replica": replica,
                        "frame": frame,
                        "time_ns": time_ns,
                        "tier": row["tier"],
                        "psi_delta": psi_delta,
                        "chi1_delta": chi1_delta,
                    }
                )
    if not selected:
        raise ValueError(f"No joint Tier A/B plus native-torsion frames in {path}")
    return selected


def cpptraj_text(
    topology: Path,
    candidates: list[dict[str, Any]],
    trajectories: dict[str, Path],
    image_anchor: str,
    output: Path,
) -> str:
    """Build cpptraj input to concatenate only eligible frames once."""
    if not candidates or not image_anchor.strip():
        raise ValueError("Candidates and an explicit image anchor are required")
    lines = [f"parm {shlex.quote(str(topology.resolve()))}"]
    for item in candidates:
        trajectory = trajectories[item["replica"]]
        frame = int(item["frame"])
        if frame < 1:
            raise ValueError("Candidate frame indices must be positive")
        source = shlex.quote(str(trajectory.resolve()))
        lines.append(f"trajin {source} {frame} {frame} 1")
    lines.extend(
        [
            f"autoimage anchor {image_anchor}",
            f"trajout {shlex.quote(str(output.resolve()))} netcdf",
            "run",
        ]
    )
    return "\n".join(lines) + "\n"


def _alignment_inputs(
    reference_cif: Path, topology_path: Path
) -> tuple[Any, Any, list[tuple[int, str]], np.ndarray]:
    """Load validated site mapping and reference fit coordinates once."""
    reference = _cif_atoms(reference_cif, "A", {auth for _, auth, _ in SITE_RESIDUES})
    fit_keys = [(auth, atom) for _, auth, _ in FIT_RESIDUES for atom in BACKBONE]
    for key in fit_keys:
        if key not in reference:
            raise ValueError(f"Missing CIF fit atom at auth residue {key[0]} {key[1]}")
    structure = parmed.load_file(str(topology_path))
    return (
        reference,
        structure,
        fit_keys,
        np.array([reference[key][1] for key in fit_keys]),
    )


def _measure_frame(
    coordinates: np.ndarray,
    reference: Any,
    structure: Any,
    fit_keys: list[tuple[int, str]],
    reference_fit: np.ndarray,
) -> dict[str, Any]:
    """Align one in-memory frame and measure the held-out native Thr233 pose."""
    if (
        coordinates.shape != (len(structure.atoms), 3)
        or not np.isfinite(coordinates).all()
    ):
        raise ValueError("Candidate coordinates do not match topology or are nonfinite")
    structure.coordinates = coordinates
    mobile = _topology_atoms(structure, SITE_RESIDUES)
    for key in fit_keys:
        if key not in mobile or reference[key][0] != mobile[key][0]:
            raise ValueError(f"Missing or element-mismatched mapped atom {key}")
    mobile_fit = np.array([mobile[key][1] for key in fit_keys])
    rotation, mobile_center, reference_center, fit_rmsd = fit_coordinates(
        mobile_fit, reference_fit
    )
    aligned = {
        key: (element, (point - mobile_center) @ rotation + reference_center)
        for key, (element, point) in mobile.items()
    }
    return {"fit_rmsd_angstrom": fit_rmsd, **state_metrics(reference, aligned)}


def _candidate_rank(
    item: dict[str, Any],
) -> tuple[float, float, float, int, float, int]:
    thr = item["native_alignment"]["thr233_held_out"]
    shifts = thr["heavy_atom_displacements_angstrom"]
    return (
        shifts["OG1"],
        thr["rmsd_angstrom"],
        shifts["N"],
        0 if item["tier"] == "A" else 1,
        item["time_ns"],
        item["frame"],
    )


def screen(
    config_path: Path, output: Path, cpptraj: Path, top_k: int
) -> dict[str, Any]:
    """Extract, align, rank and describe eligible frames without restart files."""
    if top_k < 1:
        raise ValueError("top_k must be positive")
    if output.exists():
        raise FileExistsError(output)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    topology = Path(config["topology"])
    reference_cif = Path(config["reference_cif"])
    image_anchor = str(config["image_anchor"])
    replicas = config["replicas"]
    trajectories = {item["name"]: Path(item["trajectory"]) for item in replicas}
    candidates: list[dict[str, Any]] = []
    for item in replicas:
        candidates.extend(eligible_frames(Path(item["screen_csv"]), str(item["name"])))
    for source in [
        config_path,
        topology,
        reference_cif,
        cpptraj,
        *trajectories.values(),
    ]:
        if not source.is_file():
            raise FileNotFoundError(source)
    output.mkdir(parents=True)
    compact = output / "eligible_candidates.nc"
    input_path = output / "extract_eligible.cpptraj"
    input_path.write_text(
        cpptraj_text(topology, candidates, trajectories, image_anchor, compact),
        encoding="utf-8",
    )
    amberhome = cpptraj.resolve().parent.parent
    environment = dict(os.environ)
    environment.update(AMBERHOME=str(amberhome), OMP_NUM_THREADS="1")
    library = str(amberhome / "lib")
    environment["LD_LIBRARY_PATH"] = ":".join(
        value for value in (library, environment.get("LD_LIBRARY_PATH", "")) if value
    )
    with (output / "cpptraj.log").open("w", encoding="utf-8") as log:
        subprocess.run(
            [str(cpptraj), "-i", str(input_path)],
            check=True,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    trajectory = NetCDFTraj.open_old(str(compact))
    try:
        coordinates = np.asarray(trajectory.coordinates, dtype=float)
    finally:
        trajectory.close()
    if coordinates.shape[0] != len(candidates):
        raise ValueError(
            f"Compact trajectory has {coordinates.shape[0]} frames; expected {len(candidates)}"
        )
    reference, structure, fit_keys, reference_fit = _alignment_inputs(
        reference_cif, topology
    )
    for candidate, frame_coordinates in zip(candidates, coordinates, strict=True):
        candidate["native_alignment"] = _measure_frame(
            frame_coordinates, reference, structure, fit_keys, reference_fit
        )
        thr = candidate["native_alignment"]["thr233_held_out"]
        shifts = thr["heavy_atom_displacements_angstrom"]
        candidate["accepted_native_gate"] = (
            shifts["OG1"] <= NATIVE_THR_LIMITS["og1_displacement_angstrom"]
            and shifts["N"] <= NATIVE_THR_LIMITS["n_displacement_angstrom"]
            and abs(candidate["psi_delta"])
            <= NATIVE_THR_LIMITS["absolute_psi_delta_degrees"]
            and abs(candidate["chi1_delta"])
            <= NATIVE_THR_LIMITS["absolute_chi1_delta_degrees"]
        )
    ranked = sorted(candidates, key=_candidate_rank)
    for rank, item in enumerate(ranked, start=1):
        item["rank"] = rank
    report = {
        "schema_version": 1,
        "reference_cif": str(reference_cif.resolve()),
        "alignment_fit": "backbone N/CA/C/O of Phe113, Asn186, Asp272, Asp274",
        "ranking": [
            "Thr233 OG1 displacement",
            "Thr233 heavy-atom RMSD",
            "Thr233 N displacement",
        ],
        "eligible_frame_count": len(ranked),
        "native_thr_acceptance_limits": NATIVE_THR_LIMITS,
        "accepted_count": sum(item["accepted_native_gate"] for item in ranked),
        "top_k": min(top_k, sum(item["accepted_native_gate"] for item in ranked)),
        "compact_trajectory": str(compact.resolve()),
        "sources": [
            _source_record(path)
            for path in (
                config_path,
                topology,
                reference_cif,
                cpptraj,
                *trajectories.values(),
                *(Path(item["screen_csv"]) for item in replicas),
                Path(__file__),
                input_path,
                compact,
            )
        ],
        "ranked_candidates": ranked,
        "best_candidates_for_extraction": [
            {key: item[key] for key in ("rank", "replica", "frame", "time_ns", "tier")}
            for item in [
                candidate for candidate in ranked if candidate["accepted_native_gate"]
            ][:top_k]
        ],
    }
    (output / "native_thr_screen.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (output / "ranked_candidates.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        fields = [
            "rank",
            "replica",
            "frame",
            "time_ns",
            "tier",
            "og1_displacement_angstrom",
            "thr233_heavy_atom_rmsd_angstrom",
            "n_displacement_angstrom",
            "fit_rmsd_angstrom",
            "accepted_native_gate",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in ranked:
            thr = item["native_alignment"]["thr233_held_out"]
            writer.writerow(
                {
                    "rank": item["rank"],
                    "replica": item["replica"],
                    "frame": item["frame"],
                    "time_ns": item["time_ns"],
                    "tier": item["tier"],
                    "og1_displacement_angstrom": thr[
                        "heavy_atom_displacements_angstrom"
                    ]["OG1"],
                    "thr233_heavy_atom_rmsd_angstrom": thr["rmsd_angstrom"],
                    "n_displacement_angstrom": thr["heavy_atom_displacements_angstrom"][
                        "N"
                    ],
                    "fit_rmsd_angstrom": item["native_alignment"]["fit_rmsd_angstrom"],
                    "accepted_native_gate": item["accepted_native_gate"],
                }
            )
    LOGGER.info("Aligned and ranked %d eligible WT frames", len(ranked))
    return report
