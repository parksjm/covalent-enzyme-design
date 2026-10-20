"""Export an explicit-hydrogen covalent system as a canonical-chain PDB."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def _residue_mapping() -> dict[int, tuple[str, int]]:
    mapping: dict[int, tuple[str, int]] = {}
    for first, last, chain in (
        (1, 232, "A"), (233, 464, "B"), (465, 696, "C"), (697, 928, "D"),
        (929, 1024, "E"), (1025, 1120, "F"), (1121, 1216, "G"),
        (1217, 1312, "H"),
    ):
        for amber_number in range(first, last + 1):
            mapping[amber_number] = (chain, amber_number - first + 1)
    mapping[1313] = ("I", 1)
    mapping[1314] = ("J", 1)
    return mapping


def write_inspection_pdb(topology: Any, coordinates: np.ndarray, path: Path) -> dict[int, tuple[str, int]]:
    """Write the full hydrogenated topology with canonical chain numbering.

    Args:
        topology: ParmEd Structure or AmberParm with 1314 ordered residues.
        coordinates: Cartesian coordinates in Angstrom, shape ``(atoms, 3)``.
        path: Destination PDB path.

    Returns:
        Mapping from one-based topology residue number to ``(chain, number)``.

    Raises:
        ValueError: If residue layout, coordinates, or required covalent bonds differ.
    """
    residues = list(topology.residues)
    atoms = list(topology.atoms)
    xyz = np.asarray(coordinates, dtype=float)
    if len(residues) != 1314:
        raise ValueError(f"expected 1314 residues, found {len(residues)}")
    if xyz.shape != (len(atoms), 3):
        raise ValueError(f"coordinates must have shape ({len(atoms)}, 3), got {xyz.shape}")
    if not np.isfinite(xyz).all():
        raise ValueError("coordinates contain non-finite values")
    if not all(str(residue.name).strip() for residue in residues):
        raise ValueError("all residues must have names")
    if str(residues[928].name).upper() != "THR":
        raise ValueError(f"residue 929 must be THR; found {residues[928].name}")
    if str(residues[1312].name).upper() != "LIG":
        raise ValueError(f"residue 1313 must be LIG; found {residues[1312].name}")
    if str(residues[1313].name).upper() != "WAT":
        raise ValueError(f"residue 1314 must be WAT; found {residues[1313].name}")
    mapping = _residue_mapping()
    if sorted({index for index in mapping}) != list(range(1, 1315)):
        raise ValueError("internal residue mapping is incomplete")

    atom_serial = {id(atom): index + 1 for index, atom in enumerate(atoms)}
    lig_atoms = {id(atom) for atom in residues[1312].atoms}
    og1_atoms = [atom for atom in residues[928].atoms if atom.name.strip() == "OG1"]
    if len(og1_atoms) != 1:
        raise ValueError("residue 929 must contain exactly one OG1 atom")
    bonds: list[tuple[int, int]] = []
    crosslink_found = False
    for bond in topology.bonds:
        first, second = bond.atom1, bond.atom2
        first_ligand, second_ligand = id(first) in lig_atoms, id(second) in lig_atoms
        if first_ligand and second_ligand:
            bonds.append((atom_serial[id(first)], atom_serial[id(second)]))
        elif (first is og1_atoms[0] and second_ligand) or (second is og1_atoms[0] and first_ligand):
            bonds.append((atom_serial[id(first)], atom_serial[id(second)]))
            crosslink_found = True
    if not crosslink_found:
        raise ValueError("missing Thr E1 OG1-ligand covalent bond")

    lines = ["REMARK 900 Explicit hydrogens retained; canonical chain mapping applied"]
    serial = 0
    for residue_index, residue in enumerate(residues, start=1):
        chain, residue_number = mapping[residue_index]
        for atom in residue.atoms:
            serial += 1
            x, y, z = xyz[atom_serial[id(atom)] - 1]
            name = str(atom.name).strip()
            element = str(getattr(atom, "element_name", "") or "").strip().upper()
            if not element:
                element = str(getattr(atom, "element", "")).strip().upper()
            if not element:
                raise ValueError(f"atom {name!r} has no element")
            resname = str(residue.name).strip()[:3].upper()
            lines.append(
                f"ATOM  {serial:5d} {name:>4s} {resname:3s} {chain}{residue_number:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}"
            )
        next_chain = mapping.get(residue_index + 1, (None, 0))[0]
        if next_chain != chain:
            lines.append(f"TER   {serial + 1:5d}      {str(residue.name).strip()[:3]:3s} {chain}{residue_number:4d}")
    for first, second in bonds:
        lines.append(f"CONECT{first:5d}{second:5d}")
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return mapping
