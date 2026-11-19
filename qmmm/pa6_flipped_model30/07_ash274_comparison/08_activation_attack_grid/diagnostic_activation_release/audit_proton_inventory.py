#!/usr/bin/env python3
"""Audit proton ownership and local geometry for the activation diagnostic."""

from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any
import numpy as np
import parmed

FOLDER = Path(__file__).resolve().parent
TOPOLOGY = Path(
    "/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/07_ash274_comparison/00_system/ash274_solvated.prmtop"
)
PROTONS = {
    "Thr233_N_H1": 13906,
    "Thr233_N_H2": 13907,
    "Thr233_HG1": 13917,
    "water1576_H1": 20318,
    "water1576_H2": 20319,
    "ASH274_HD2": 14580,
    "substrate_leaving_HN": 19941,
}
ATOMS = {
    "Thr233_N": 13905,
    "Thr233_N_H1": 13906,
    "Thr233_N_H2": 13907,
    "Thr233_OG1": 13916,
    "Thr233_HG1": 13917,
    "substrate_C": 19920,
    "substrate_O": 19921,
    "substrate_leaving_N": 19933,
    "ASH274_OD2": 14579,
    "ASH274_HD2": 14580,
    "water1576_O": 20317,
    "water1576_H1": 20318,
    "water1576_H2": 20319,
}


def digest(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_geometry(restart: Path) -> tuple[Any, np.ndarray]:
    """Load an Amber restart with its topology and return coordinates."""
    system = parmed.load_file(str(TOPOLOGY), xyz=str(restart))
    return system, np.asarray(system.coordinates, dtype=float)


def snapshot(restart: Path) -> dict[str, Any]:
    """Measure nearest-heavy ownership and required local distances."""
    system, xyz = load_geometry(restart)
    heavy = [atom.idx for atom in system.atoms if atom.atomic_number > 1]
    inventory: dict[str, Any] = {}
    for label, atom_id in PROTONS.items():
        h_idx = atom_id - 1
        distances = np.linalg.norm(xyz[heavy] - xyz[h_idx], axis=1)
        nearest = int(np.argmin(distances))
        owner = system.atoms[heavy[nearest]]
        inventory[label] = {
            "nearest_heavy_atom_id": owner.idx + 1,
            "nearest_heavy_atom_name": owner.name,
            "nearest_heavy_residue": owner.residue.name,
            "nearest_heavy_residue_number": owner.residue.number,
            "distance_A": float(distances[nearest]),
        }

    def dist(first: str, second: str) -> float:
        return float(np.linalg.norm(xyz[ATOMS[first] - 1] - xyz[ATOMS[second] - 1]))

    geometry = {
        "q_N_A": dist("Thr233_OG1", "Thr233_HG1") - dist("Thr233_N", "Thr233_HG1"),
        "attack_OG1_substrate_C_A": dist("Thr233_OG1", "substrate_C"),
        "substrate_C_O_A": dist("substrate_C", "substrate_O"),
        "substrate_C_leaving_N_A": dist("substrate_C", "substrate_leaving_N"),
        "Thr233_N_H1_A": dist("Thr233_N", "Thr233_N_H1"),
        "Thr233_N_H2_A": dist("Thr233_N", "Thr233_N_H2"),
        "Thr233_N_HG1_A": dist("Thr233_N", "Thr233_HG1"),
        "Thr233_OG1_HG1_A": dist("Thr233_OG1", "Thr233_HG1"),
        "ASH274_OD2_HD2_A": dist("ASH274_OD2", "ASH274_HD2"),
        "water1576_O_H1_A": dist("water1576_O", "water1576_H1"),
        "water1576_O_H2_A": dist("water1576_O", "water1576_H2"),
    }
    return {
        "restart_sha256": digest(restart),
        "proton_nearest_heavy": inventory,
        "geometry": geometry,
    }


def main() -> None:
    """Write a hashed initial/endpoint proton-inventory report."""
    source = FOLDER / "starting.ncrst"
    endpoint = FOLDER / "optimized.ncrst"
    if not endpoint.is_file():
        raise FileNotFoundError(endpoint)
    report = {
        "schema_version": 1,
        "classification": "diagnostic only; no validated TI1 or physical-intermediate claim",
        "topology_sha256": digest(TOPOLOGY),
        "initial": snapshot(source),
        "endpoint": snapshot(endpoint),
    }
    (FOLDER / "proton_inventory.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
