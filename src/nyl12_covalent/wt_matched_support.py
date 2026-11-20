"""Coordinate-preserving helpers for matched WT protonation preparation."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np
from MDAnalysis.lib.distances import distance_array

BETA = (929, 1025, 1121, 1217)
ASH = (970, 1066, 1162, 1258)
ENDS = (232, 464, 696, 928, 1024, 1120, 1216, 1312)
SOLVENT = {"WAT", "Na+", "Cl-"}


def digest(path: Path) -> str:
    """Hash an existing file without changing it."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distant_sodium(parent: Any) -> list[int]:
    """Select four sodium residues farthest from solute under periodic MIC."""
    sodium = [a for a in parent.atoms if a.residue.name == "Na+"]
    heavy = [
        a.idx
        for a in parent.atoms
        if a.residue.name not in SOLVENT and a.atomic_number > 1
    ]
    distances = distance_array(
        parent.coordinates[[a.idx for a in sodium]],
        parent.coordinates[heavy],
        box=parent.box,
    )
    scores = distances.min(axis=1)
    ranked = sorted(
        zip(scores, sodium, strict=True), key=lambda x: (-float(x[0]), x[1].residue.idx)
    )
    if len(ranked) < 4:
        raise ValueError("Need at least four sodium ions")
    return [atom.residue.idx + 1 for _, atom in ranked[:4]]


def write_full_pdb(parent: Any, path: Path, removed: list[int], ash: bool) -> list[int]:
    """Write all retained residues; return output-to-parent residue IDs."""
    coordinates = parent.coordinates
    mapping = []
    lines = []
    serial = 0
    for residue in parent.residues:
        original = residue.idx + 1
        if original in removed:
            continue
        mapping.append(original)
        name = "ASH" if ash and original in ASH else residue.name
        for atom in residue.atoms:
            serial += 1
            x, y, z = coordinates[atom.idx]
            lines.append(
                f"ATOM  {serial % 100000:5d} {atom.name:>4} {name:>3} "
                f" {len(mapping) % 10000:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
            )
        if original in ENDS:
            lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return mapping


def bond_keys(
    structure: Any, residue_map: list[int]
) -> set[tuple[tuple[int, str], ...]]:
    """Represent bonds by stable parent residue IDs and atom names."""
    return {
        tuple(
            sorted(
                (
                    (residue_map[b.atom1.residue.idx], b.atom1.name),
                    (residue_map[b.atom2.residue.idx], b.atom2.name),
                )
            )
        )
        for b in structure.bonds
    }


def restore_and_check(
    parent: Any, product: Any, mapping: list[int], charged: bool, ash: bool
) -> dict:
    """Restore shared coordinates and validate graph, charge and new protons."""
    if len(product.residues) != len(mapping):
        raise ValueError("LEaP changed residue count/order")
    expected_new = {(r, "H3") for r in BETA} if charged else set()
    if ash:
        expected_new |= {(r, "HD2") for r in ASH}
    parent_xyz = parent.coordinates
    product_xyz = product.coordinates
    found_new = set()
    shared = set()
    for residue, original in zip(product.residues, mapping, strict=True):
        old = parent.residues[original - 1]
        expected_name = "ASH" if ash and original in ASH else old.name
        if residue.name != expected_name:
            raise ValueError(f"Residue mapping changed {original}")
        old_atoms = {a.name: a for a in old.atoms}
        for atom in residue.atoms:
            key = (original, atom.name)
            if atom.name in old_atoms:
                product_xyz[atom.idx] = parent_xyz[old_atoms[atom.name].idx]
                shared.add(key)
            else:
                found_new.add(key)
    if found_new != expected_new:
        raise ValueError(f"Unexpected added H inventory {found_new}")
    wanted = {(r, a.name) for r in mapping for a in parent.residues[r - 1].atoms}
    if shared != wanted:
        raise ValueError("A shared parent atom was lost")
    old_graph = bond_keys(parent, list(range(1, len(parent.residues) + 1)))
    new_graph = bond_keys(product, mapping)
    retained = {b for b in old_graph if all(k in shared for k in b)}
    actual = {b for b in new_graph if all(k in shared for k in b)}
    if actual != retained:
        raise ValueError(f"Shared bond graph changed: {actual ^ retained}")
    expected_bonds = {
        tuple(sorted((key, (key[0], "N" if key[1] == "H3" else "OD2"))))
        for key in expected_new
    }
    if new_graph - actual != expected_bonds:
        raise ValueError("New proton connectivity mismatch")
    product.coordinates = product_xyz
    coords = product_xyz.copy()
    if not np.isfinite(coords).all():
        raise ValueError("Nonfinite coordinates")
    new_ids = [
        a.idx for a in product.atoms if (mapping[a.residue.idx], a.name) in expected_new
    ]
    if new_ids:
        ds = distance_array(coords[new_ids], coords, box=parent.box)
        for row, idx in enumerate(new_ids):
            ds[row, idx] = np.inf
        if ds.min() < 0.45:
            raise ValueError("New proton superposition")
    dry_charge = sum(a.charge for a in product.atoms if a.residue.name not in SOLVENT)
    target = -4 if charged or ash else -8
    if abs(dry_charge - target) > 1e-5:
        raise ValueError(f"Incorrect dry charge {dry_charge}")
    for r in BETA:
        atoms = {a.name: a for a in product.residues[mapping.index(r)].atoms}
        if ("H3" in atoms) != charged or not {"H1", "H2", "HG1"} <= atoms.keys():
            raise ValueError("Incorrect Thr hydrogen inventory")
        if abs(sum(a.charge for a in atoms.values()) - int(charged)) > 1e-5:
            raise ValueError("Incorrect whole Thr charge")
    return {
        "dry_charge": dry_charge,
        "expected_dry_charge": target,
        "shared_atoms_exact": len(shared),
        "new_hydrogens": sorted(expected_new),
        "shared_bond_graph_exact": True,
        "new_proton_min_distance_A": float(ds.min()) if new_ids else None,
    }
