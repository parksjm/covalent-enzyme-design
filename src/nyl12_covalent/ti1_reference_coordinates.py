"""Construct constrained TI1 surrogate seeds from a named protein reference."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np
from parmed.formats import CIFFile
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

TORSION_LABELS = {
    "phi": ("ACE:C", "THR:N", "THR:CA", "THR:C"),
    "psi": ("THR:N", "THR:CA", "THR:C", "NME:N"),
    "chi1": ("THR:N", "THR:CA", "THR:CB", "THR:OG1"),
}


def _dihedral(points: list[np.ndarray]) -> float:
    conformer = Chem.Conformer(4)
    for index, point in enumerate(points):
        conformer.SetAtomPosition(index, tuple(float(x) for x in point))
    value = float(rdMolTransforms.GetDihedralDeg(conformer, 0, 1, 2, 3))
    if not np.isfinite(value):
        raise ValueError("Degenerate source dihedral")
    return value


def crystal_backbone_reference(
    path: Path, chain: str, residue_number: int
) -> dict[str, Any]:
    """Read explicitly selected Thr torsions from a single-model mmCIF.

    Args:
        path: Unmodified reference mmCIF.
        chain: Exact author-chain identifier parsed by ParmEd.
        residue_number: Exact author residue number of catalytic Thr.

    Returns:
        Source identity, selected coordinates, psi/chi1, and phi only when
        the preceding carbonyl is present at peptide-bond distance.

    Raises:
        ValueError: If selection, atom inventory, or geometry is ambiguous.
    """
    structure = CIFFile.parse(str(path), skip_bonds=True)
    if structure.get_coordinates().shape[0] != 1:
        raise ValueError("Select a single-model reference before preparation")

    def residue_atoms(number: int) -> dict[str, np.ndarray]:
        residues = [
            r for r in structure.residues if r.chain == chain and r.number == number
        ]
        if len(residues) != 1:
            raise ValueError(f"Expected one {chain}:{number}, found {len(residues)}")
        residue = residues[0]
        if number == residue_number and residue.name != "THR":
            raise ValueError(f"{chain}:{number} is {residue.name}, not THR")
        names = [a.name for a in residue.atoms]
        if len(names) != len(set(names)):
            raise ValueError("Ambiguous duplicate source atom names")
        if any(a.other_locations for a in residue.atoms):
            raise ValueError("Resolve alternate conformations explicitly")
        return {
            a.name: np.asarray((a.xx, a.xy, a.xz), dtype=float) for a in residue.atoms
        }

    current = residue_atoms(residue_number)
    following = residue_atoms(residue_number + 1)
    required = {"N", "CA", "C", "O", "CB", "OG1", "CG2"}
    if required - current.keys() or "N" not in following:
        raise ValueError("Source lacks required Thr/next-residue atoms")
    if not all(np.isfinite(v).all() for v in current.values()):
        raise ValueError("Source contains nonfinite coordinates")
    next_distance = float(np.linalg.norm(current["C"] - following["N"]))
    if not 1.1 <= next_distance <= 1.7:
        raise ValueError("Selected next N is not at peptide-bond distance")
    psi = _dihedral([current[x] for x in ("N", "CA", "C")] + [following["N"]])
    chi1 = _dihedral([current[x] for x in ("N", "CA", "CB", "OG1")])
    preceding = [
        r
        for r in structure.residues
        if r.chain == chain and r.number == residue_number - 1
    ]
    phi, previous_distance = None, None
    if preceding:
        previous = residue_atoms(residue_number - 1)
        if "C" not in previous:
            raise ValueError("Preceding residue lacks carbonyl C")
        previous_distance = float(np.linalg.norm(previous["C"] - current["N"]))
        if 1.1 <= previous_distance <= 1.7:
            phi = _dihedral([previous["C"]] + [current[x] for x in ("N", "CA", "C")])
    return {
        "source": str(path.absolute()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "chain": chain,
        "residue_number": residue_number,
        "phi_degrees": phi,
        "phi_scope": "preceding-peptide reference" if phi is not None else "undefined",
        "psi_degrees": psi,
        "chi1_degrees": chi1,
        "previous_C_to_N_A": previous_distance,
        "C_to_next_N_A": next_distance,
        "thr_coordinates_A": {
            name: current[name].tolist() for name in sorted(required)
        },
    }


def embed_bonded_reference(
    molecule: Chem.Mol,
    atom_indices: dict[str, int],
    targets: dict[str, float],
    seed: int,
) -> Chem.Mol:
    """Embed one deterministic seed and set the three declared torsions.

    Args:
        molecule: Explicit-hydrogen, stereochemically defined surrogate graph.
        atom_indices: Unique residue-qualified atom indexes.
        targets: Explicit phi, psi and chi1 angles in degrees.
        seed: Reproducible ETKDG random seed.

    Returns:
        A copy with one conformer; no force-field minimization is performed.

    Raises:
        ValueError: If inputs, embedding, torsions or chirality fail.
    """
    if set(targets) != set(TORSION_LABELS):
        raise ValueError("Exactly phi, psi and chi1 targets are required")
    if not all(np.isfinite(v) for v in targets.values()):
        raise ValueError("Nonfinite torsion target")
    result = Chem.Mol(molecule)
    result.RemoveAllConformers()
    settings = AllChem.ETKDGv3()
    settings.randomSeed = seed
    settings.enforceChirality = True
    if AllChem.EmbedMolecule(result, settings) != 0:
        raise ValueError("ETKDG seed construction failed")
    expected = {
        a.GetIdx(): a.GetProp("_CIPCode")
        for a in result.GetAtoms()
        if a.HasProp("_CIPCode")
    }
    conformer = result.GetConformer()
    for name, labels in TORSION_LABELS.items():
        indices = [atom_indices[label] for label in labels]
        rdMolTransforms.SetDihedralDeg(conformer, *indices, targets[name])
    Chem.AssignStereochemistryFrom3D(result, replaceExistingTags=True)
    for index, target in expected.items():
        atom = result.GetAtomWithIdx(index)
        if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != target:
            raise ValueError("Seed construction changed stereochemistry")
    return result
