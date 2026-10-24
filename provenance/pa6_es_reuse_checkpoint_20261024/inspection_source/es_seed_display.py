"""Mapped ES seed presentation and rigid alignment helpers."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
from rdkit.Numerics import rdAlignment


def _display(number: int, water_numbers: dict[int, int]) -> tuple[str, int]:
    if number <= 928:
        copy, local = divmod(number - 1, 232)
        return "ABCD"[copy], local + 1
    if number <= 1312:
        copy, local = divmod(number - 929, 96)
        return "EFGH"[copy], local + 233
    if number == 1313:
        return "I", 1
    return "W", water_numbers[number]


def _write_pdb(parm: Any, residues: set[int], waters: list[int], path: Path) -> dict:
    xyz = parm.coordinates.copy()
    water_numbers = {number: i + 1 for i, number in enumerate(waters)}
    atoms = [atom for atom in parm.atoms if atom.residue.idx + 1 in residues]
    serials = {atom.idx: i + 1 for i, atom in enumerate(atoms)}
    records, lines = [], ["REMARK ES SEED; NOT AN OPTIMIZED TI OR TS\n"]
    for atom in atoms:
        number = atom.residue.idx + 1
        chain, display = _display(number, water_numbers)
        x, y, z = xyz[atom.idx]
        if not np.isfinite([x, y, z]).all() or max(abs(x), abs(y), abs(z)) >= 1000:
            raise ValueError("Coordinates cannot be represented safely in PDB")
        serial = serials[atom.idx]
        element = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}[atom.atomic_number]
        record = "ATOM" if number <= 1312 else "HETATM"
        lines.append(
            f"{record:<6}{serial:5d} {atom.name:^4s} {atom.residue.name:>3s} "
            f"{chain}{display:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        )
        records.append(
            {
                "serial": serial,
                "original_atom_index": atom.idx + 1,
                "original_residue_index": number,
                "atom_name": atom.name,
                "display_chain": chain,
                "display_residue": display,
            }
        )
    for bond in parm.bonds:
        if bond.atom1.idx in serials and bond.atom2.idx in serials:
            lines.append(
                f"CONECT{serials[bond.atom1.idx]:5d}{serials[bond.atom2.idx]:5d}\n"
            )
    path.write_text("".join(lines) + "END\n")
    return {"path": str(path), "atoms": records}


def _align(parm: Any, reference: np.ndarray | None) -> tuple[np.ndarray, dict]:
    indices = []
    for residue in parm.residues[:1312]:
        atoms = [atom.idx for atom in residue.atoms if atom.name == "CA"]
        if len(atoms) != 1:
            raise ValueError("Expected one CA in every protein residue")
        indices.extend(atoms)
    if len(indices) != 1312:
        raise ValueError("Expected all1312 protein CA atoms")
    xyz = parm.coordinates.copy()
    probe = xyz[indices].copy()
    if reference is None:
        reference = probe.copy()
    _, transform = rdAlignment.GetAlignmentTransform(reference, probe)
    xyz = xyz @ transform[:3, :3].T + transform[:3, 3]
    parm.coordinates = xyz
    rmsd = float(np.sqrt(np.mean(np.sum((xyz[indices] - reference) ** 2, axis=1))))
    return reference, {
        "transform": transform.tolist(),
        "protein_ca_rmsd_A": rmsd,
        "protein_ca_atom_indices": [index + 1 for index in indices],
    }


def _viewer(output: Path, names: list[str]) -> None:
    commands = [
        "from pathlib import Path",
        "from pymol import cmd, util",
        "base = Path(__file__).resolve().parent",
        "cmd.load(str(base / 'fullprotein_substrate.pdb'), 'context')",
    ]
    commands += [f"cmd.load(str(base / '{name}.pdb'), '{name}')" for name in names]
    commands += [
        "cmd.hide('everything', 'all')",
        "cmd.show('cartoon', 'context')",
        "cmd.color('gray80', 'context')",
        "cmd.set('cartoon_transparency', 0.65)",
        "cmd.show('sticks', 'not context')",
        "cmd.show('spheres', 'chain W and elem H')",
        "cmd.set('sphere_scale', 0.18, 'chain W and elem H')",
        "cmd.show('sticks', '(resi 113+186 and chain A) and not context')",
        "cmd.set('valence', 0)",
        "util.cnc()",
    ]
    commands += [f"cmd.disable('{name}')" for name in names[1:]]
    commands += [f"cmd.zoom('{names[0]}', 5)"]
    (output / "view_es_seeds.py").write_text("\n".join(commands) + "\n")
    (output / "view_es_seeds.pml").write_text(
        "# First cd to this inspection folder in PyMOL.\nrun view_es_seeds.py\n"
    )
