"""Audit and freeze the corrected ASH274 movable-shell definition."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from parmed.amber.mask import AmberMask
from scipy.spatial import cKDTree

CASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
TOPOLOGY = CASE.parent.parent / "00_system/ash274_solvated.prmtop"
OLD_PREPARATION = CASE / "preparation.json"
RESTARTS = {
    "original_starting": CASE / "starting.ncrst",
    "prior_optimized": CASE / "optimized.ncrst",
}
RADIUS_A = 12.0


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def residue_mask(residues: set[int]) -> str:
    """Compress one-based residue numbers into an Amber mask."""
    values = sorted(residues)
    ranges: list[str] = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ":" + ",".join(ranges)


def main() -> None:
    """Write a provenance-rich, non-overwritable shell audit."""
    output = OUT / "shell_audit.json"
    if output.exists():
        raise FileExistsError(output)
    preparation = json.loads(OLD_PREPARATION.read_text())
    topology = parmed.load_file(str(TOPOLOGY))
    qm_mask = preparation["qm_region"]["mask"]
    old_mask = preparation["movable_region"]["mask"]
    qm_indices = [
        i
        for i, selected in enumerate(AmberMask(topology, qm_mask).Selection())
        if selected
    ]
    old_indices = [
        i
        for i, selected in enumerate(AmberMask(topology, old_mask).Selection())
        if selected
    ]
    if not qm_indices or not old_indices:
        raise ValueError("QM or inherited active mask selected no atoms")
    old_residues = {int(topology.atoms[i].residue.idx) + 1 for i in old_indices}
    if len(old_indices) != preparation["movable_region"]["atom_count"]:
        raise ValueError("Inherited active mask does not reproduce recorded atom count")

    per_restart: dict[str, Any] = {}
    shell_residues = set(old_residues)
    for label, restart in RESTARTS.items():
        structure = parmed.load_file(str(TOPOLOGY), xyz=str(restart))
        if len(structure.atoms) != len(topology.atoms):
            raise ValueError(f"Atom count mismatch in {restart}")
        coordinates = np.asarray(structure.coordinates, dtype=float)
        qm_tree = cKDTree(coordinates[qm_indices])
        distances, _ = qm_tree.query(coordinates, k=1)
        near_atoms = set(np.flatnonzero(distances <= RADIUS_A).tolist())
        near_residues = {int(topology.atoms[i].residue.idx) + 1 for i in near_atoms}
        omitted_residues = near_residues - old_residues
        shell_residues.update(near_residues)
        water1576 = [atom.idx for atom in topology.residues[1575].atoms]
        per_restart[label] = {
            "path": str(restart.resolve()),
            "sha256": sha256(restart),
            "atoms_within_12A": len(near_atoms),
            "residues_within_12A": len(near_residues),
            "atoms_within_12A_outside_old_mask": len(near_atoms - set(old_indices)),
            "residues_within_12A_outside_old_mask": len(omitted_residues),
            "omitted_residue_numbers": sorted(omitted_residues),
            "water_1576_OG1_distance_A": float(
                np.linalg.norm(
                    coordinates[water1576[0]]
                    - coordinates[
                        next(
                            atom.idx
                            for atom in topology.residues[928].atoms
                            if atom.name == "OG1"
                        )
                    ]
                )
            ),
            "water_1576_minimum_distance_to_qm_A": float(np.min(distances[water1576])),
        }
    new_mask = residue_mask(shell_residues)
    new_indices = [
        i
        for i, selected in enumerate(AmberMask(topology, new_mask).Selection())
        if selected
    ]
    if set(qm_indices) - set(new_indices):
        raise ValueError("Corrected shell excludes QM atom(s)")
    if not set(old_indices) <= set(new_indices):
        raise ValueError("Corrected shell does not contain inherited selected atoms")
    if 1576 not in shell_residues:
        raise ValueError("Water residue 1576 is not movable in corrected shell")
    report = {
        "purpose": "Correct inherited ASH274 E-S DFTB3/MM active shell; no QM region or chemistry change.",
        "method": "cKDTree nearest Euclidean atom distance, evaluated independently on both recorded restart coordinate sets; complete residues selected.",
        "radius_angstrom": RADIUS_A,
        "topology": {
            "path": str(TOPOLOGY.resolve()),
            "sha256": sha256(TOPOLOGY),
            "atoms": len(topology.atoms),
            "residues": len(topology.residues),
        },
        "old_preparation": {
            "path": str(OLD_PREPARATION.resolve()),
            "sha256": sha256(OLD_PREPARATION),
        },
        "qm": {
            "mask": qm_mask,
            "charge": preparation["qm_region"]["charge"],
            "atom_count": len(qm_indices),
        },
        "old_active": {
            "mask": old_mask,
            "atom_count": len(old_indices),
            "residue_count": len(old_residues),
            "residues": sorted(old_residues),
        },
        "restart_audits": per_restart,
        "new_active": {
            "mask": new_mask,
            "atom_count": len(new_indices),
            "residue_count": len(shell_residues),
            "added_residues": sorted(shell_residues - old_residues),
            "water_1576_included": True,
        },
        "preserved_protocol": {
            "charge": preparation["qm_region"]["charge"],
            "proton_inventory": preparation["microstate"],
            "periodic_embedding": preparation["periodic_embedding"],
            "optimization": preparation["optimization"],
        },
        "source_hashes": {
            str(path.resolve()): sha256(path)
            for path in [Path(__file__), TOPOLOGY, OLD_PREPARATION, *RESTARTS.values()]
        },
    }
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
