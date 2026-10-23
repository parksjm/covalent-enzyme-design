"""Extract source-mapped neutral amide donors with direction-aware C-H caps."""

from __future__ import annotations

import numpy as np
import parmed as pmd
from rdkit import Chem


def _donor(
    model: pmd.Structure,
    prefix: str,
    selected: list[tuple[int, str]],
    cuts: list[tuple[int, str, str, str]],
) -> tuple[Chem.Mol, list, list]:
    rw, points, records, source_indices = Chem.RWMol(), [], [], {}
    for residue_number, name in selected:
        residue = model.residues[residue_number - 1]
        matches = [a for a in residue.atoms if a.name == name]
        if len(matches) != 1:
            raise ValueError(f"Missing/duplicate donor atom {residue_number}:{name}")
        source = matches[0]
        atom = Chem.Atom(source.atomic_number)
        atom.SetNoImplicit(True)
        origin = f"{prefix}:{residue_number}:{name}"
        atom.SetProp("origin_id", origin)
        source_indices[source.idx] = rw.AddAtom(atom)
        points.append(model.coordinates[source.idx].copy())
        records.append(
            {
                "origin_id": origin,
                "source_atom_index_1based": source.idx + 1,
                "source_residue_1based": residue_number,
                "source_residue_name": residue.name,
                "source_atom_name": name,
                "source_coordinates": "prepared.rst7",
                "constructed_cap_H": False,
            }
        )
    for bond in model.bonds:
        left, right = bond.atom1, bond.atom2
        if left.idx in source_indices and right.idx in source_indices:
            double = {left.name, right.name} in ({"C", "O"}, {"CG", "OD1"})
            rw.AddBond(
                source_indices[left.idx],
                source_indices[right.idx],
                Chem.BondType.DOUBLE if double else Chem.BondType.SINGLE,
            )
    for residue_number, center_name, removed_name, cap_name in cuts:
        atoms = {a.name: a for a in model.residues[residue_number - 1].atoms}
        center, removed = atoms[center_name], atoms[removed_name]
        if removed not in center.bond_partners:
            raise ValueError("Declared cap direction is not a source bond")
        direction = model.coordinates[removed.idx] - model.coordinates[center.idx]
        if np.linalg.norm(direction) < 1e-8:
            raise ValueError("Degenerate cap direction")
        atom = Chem.Atom(1)
        atom.SetNoImplicit(True)
        origin = f"{prefix}:{residue_number}:{cap_name}"
        atom.SetProp("origin_id", origin)
        index = rw.AddAtom(atom)
        rw.AddBond(source_indices[center.idx], index, Chem.BondType.SINGLE)
        points.append(
            model.coordinates[center.idx] + 1.09 * direction / np.linalg.norm(direction)
        )
        records.append(
            {
                "origin_id": origin,
                "constructed_cap_H": True,
                "source_coordinates": "prepared.rst7",
                "CH_distance_A": 1.09,
                "cut_source_atom_indices_1based": [center.idx + 1, removed.idx + 1],
                "cut": [
                    f"{residue_number}:{center_name}",
                    f"{residue_number}:{removed_name}",
                ],
            }
        )
    molecule = rw.GetMol()
    Chem.SanitizeMol(molecule)
    return molecule, points, records
