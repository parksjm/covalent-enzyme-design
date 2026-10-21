"""Geometry and stage-specific energy helpers for Gaussian conformer audits."""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path

import numpy as np
from rdkit import Chem


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_stage_energies(text: str) -> dict[str, float]:
    """Return optimized DFT and HF ESP SCF energies from their Link1 jobs."""
    blocks = re.split(r"Normal termination of Gaussian 16", text)
    optimized = [block for block in blocks[:-1] if "Optimization completed." in block]
    if not optimized:
        raise ValueError("Expected a normally terminated optimization job")
    number = r"([-+]?\d+(?:\.\d*)?(?:[DEde][-+]?\d+)?)"
    dft_matches = re.findall(
        r"SCF Done:\s+E\(RwB97XD\)\s*=\s*" + number,
        optimized[-1].split("Optimization completed.", 1)[0],
    )
    after_optimization = blocks[blocks.index(optimized[-1]) + 1 : -1]
    hf_blocks = [block for block in after_optimization if "E(RHF)" in block]
    if len(hf_blocks) != 1:
        raise ValueError(
            "Expected one normally terminated HF ESP job after optimization"
        )
    hf_matches = re.findall(r"SCF Done:\s+E\(RHF\)\s*=\s*" + number, hf_blocks[0])
    if not dft_matches or len(hf_matches) != 1:
        raise ValueError(
            "Could not uniquely identify optimized DFT and HF ESP energies"
        )
    values = {
        "optimized_wb97xd_Eh": float(
            dft_matches[-1].replace("D", "E").replace("d", "e")
        ),
        "hf_esp_Eh": float(hf_matches[0].replace("D", "E").replace("d", "e")),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("Gaussian energy is nonfinite")
    return values


def heavy_atom_kabsch_rmsd(
    first: np.ndarray, second: np.ndarray, elements: list[str]
) -> float:
    """Compute same-order heavy-atom RMSD after proper-rotation Kabsch alignment."""
    indices = [index for index, element in enumerate(elements) if element != "H"]
    left, right = first[indices], second[indices]
    left_centered = left - left.mean(axis=0)
    right_centered = right - right.mean(axis=0)
    u, _, vt = np.linalg.svd(left_centered.T @ right_centered)
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    aligned = left_centered @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - right_centered) ** 2, axis=1))))


def dihedral_degrees(points: np.ndarray) -> float:
    """Calculate a signed four-point torsion in degrees."""
    b0 = -(points[1] - points[0])
    b1 = points[2] - points[1]
    b2 = points[3] - points[2]
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def named_psi_indices(atoms: list[dict[str, object]]) -> list[int]:
    """Resolve mapped atom indices for N-CA-C-NME_N psi."""
    wanted = (("THR", "N"), ("THR", "CA"), ("THR", "C"), ("NME", "N"))
    result = []
    for residue, name in wanted:
        matches = [
            int(atom["index"])
            for atom in atoms
            if atom.get("residue") == residue and atom.get("name") == name
        ]
        if len(matches) != 1:
            raise ValueError(f"Atom map does not uniquely identify {residue}:{name}")
        result.append(matches[0])
    return result


def write_optimized_sdf(molecule_path: Path, xyz: np.ndarray, output: Path) -> None:
    """Write original molecular graph with audited final coordinates."""
    molecules = list(Chem.SDMolSupplier(str(molecule_path), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError(
            f"Expected one explicit-hydrogen SDF molecule: {molecule_path}"
        )
    molecule = molecules[0]
    if molecule.GetNumAtoms() != len(xyz):
        raise ValueError("Optimized coordinate count differs from SDF atom inventory")
    molecule.RemoveAllConformers()
    conformer = Chem.Conformer(len(xyz))
    conformer.Set3D(True)
    for index, point in enumerate(xyz):
        conformer.SetAtomPosition(index, tuple(float(value) for value in point))
    molecule.AddConformer(conformer)
    molecule.SetProp("coordinate_source", "final Gaussian optimized orientation")
    molecule.SetProp("coordinate_units", "angstrom")
    molecule.SetProp("source_molecule_sha256", sha256(molecule_path))
    with Chem.SDWriter(str(output)) as writer:
        writer.write(molecule)
