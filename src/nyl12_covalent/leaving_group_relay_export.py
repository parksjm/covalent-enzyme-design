"""Native geometry reporting and exact source-coordinate relay exports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from parmed.amber import NetCDFRestart
from scipy.io import netcdf_file

from nyl12_covalent.catalytic_alignment import SITE_RESIDUES, _cif_atoms
from nyl12_covalent.leaving_group_relay_mapping import resolve_atom
from nyl12_covalent.nalpha_water_relay import cell_matrix, hbond, native_geometry


def digest(path: Path) -> str:
    """Return streaming source SHA-256."""
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def source_frame(row: dict) -> tuple[np.ndarray, np.ndarray]:
    """Load the exact indexed source frame and verify its recorded time."""
    with netcdf_file(row["trajectory"], mmap=True) as source:
        index = row["local_frame"] - 1
        if abs(float(source.variables["time"][index]) - row["time_ps"]) > 1e-3:
            raise ValueError("Source frame time mismatch")
        xyz = np.array(source.variables["coordinates"][index], dtype=float)
        box = np.concatenate(
            [
                np.array(source.variables[k][index])
                for k in ("cell_lengths", "cell_angles")
            ]
        )
    return xyz, box


def evaluate_native(rows: list[dict], config: dict) -> list[dict]:
    """Evaluate only supplied frames against the existing native Thr gates."""
    expected = {
        "og1_shift_A": 0.75,
        "nalpha_shift_A": 1.25,
        "abs_psi_delta_deg": 35,
        "abs_chi1_delta_deg": 20,
    }
    if config["native_thresholds"] != expected:
        raise ValueError("Configured native thresholds differ from established gates")
    reference = _cif_atoms(
        Path(config["reference_cif"]), "A", {auth for _, auth, _ in SITE_RESIDUES}
    )
    cases = {case["name"]: case for case in config["cases"]}
    topologies: dict[str, Any] = {}
    result = []
    for row in rows:
        case = cases[row["case"]]
        if row["case"] not in topologies:
            topologies[row["case"]] = parmed.load_file(case["topology"])
        topology = topologies[row["case"]]
        xyz, box = source_frame(row)
        nalpha = resolve_atom(topology, case["geometry_atoms"]["nalpha"])
        metrics = native_geometry(
            xyz, topology, reference, cell_matrix(box[:3], box[3:]), nalpha
        )
        activation = {}
        if row.get("nalpha_NH2_start", False):
            og1 = resolve_atom(topology, case["geometry_atoms"]["og1"])
            hydrogens = [
                a.idx for a in topology.atoms[og1].bond_partners if a.atomic_number == 1
            ]
            if not hydrogens:
                raise ValueError(
                    "OG1 activation assessment requires its explicit owned hydrogen"
                )
            legs = [
                hbond(xyz, og1, h, nalpha, cell_matrix(box[:3], box[3:]))
                for h in hydrogens
            ]
            leg = min(legs, key=lambda edge: (not edge["pass"], edge["ha_A"]))
            activation = {
                "og1_to_nalpha_activation_leg": leg,
                "activation_status": "direct_activation_observed"
                if leg["pass"]
                else "activation_unresolved",
            }
        failures = [
            name
            for name, passed in {
                "og1_shift": metrics["native_og1_shift_A"] <= expected["og1_shift_A"],
                "nalpha_shift": metrics["native_nalpha_shift_A"]
                <= expected["nalpha_shift_A"],
                "psi": abs(metrics["psi_delta_deg"]) <= expected["abs_psi_delta_deg"],
                "chi1": abs(metrics["chi1_delta_deg"])
                <= expected["abs_chi1_delta_deg"],
            }.items()
            if not passed
        ]
        result.append({**row, **metrics, **activation, "native_failed_gates": failures})
    return result


def export_frame(row: dict, case: dict, output: Path, label: str) -> dict:
    """Export and exactly verify a full-system source frame, with atom provenance."""
    topology = parmed.load_file(case["topology"])
    xyz, box = source_frame(row)
    if xyz.shape != (len(topology.atoms), 3):
        raise ValueError("Topology/source atom count mismatch")
    path = output / f"{label}.ncrst"
    restart = NetCDFRestart.open_new(str(path), len(xyz), True, False, title=label)
    try:
        restart.coordinates, restart.box = xyz, box
    finally:
        restart.close()
    verify = NetCDFRestart.open_old(str(path))
    try:
        if not np.array_equal(xyz, np.asarray(verify.coordinates).reshape(-1, 3)):
            raise ValueError("Export changed source coordinates")
        if not np.array_equal(box, np.asarray(verify.box).reshape(-1)):
            raise ValueError("Export changed source box")
    finally:
        verify.close()
    water_ids = sorted(
        int(label.split(":")[1]) for label in row["path"] if label.startswith("water:")
    )
    record = {
        **row,
        "restart": str(path.resolve()),
        "restart_sha256": digest(path),
        "topology": case["topology"],
        "topology_sha256": digest(Path(case["topology"])),
        "coordinates_unmodified": True,
        "box_unmodified": True,
        "coordinates_float64_sha256": hashlib.sha256(
            xyz.astype("<f8").tobytes()
        ).hexdigest(),
        "water_residue_indices1": water_ids,
        "water_atoms_index1": {
            str(i): [a.idx + 1 for a in topology.residues[i - 1].atoms]
            for i in water_ids
        },
        "geometry_atoms_index1": {
            name: resolve_atom(topology, descriptor) + 1
            for name, descriptor in case["geometry_atoms"].items()
        },
        "velocities_included": False,
    }
    (output / f"{label}.json").write_text(json.dumps(record, indent=2) + "\n")
    return record
