"""Generate traceable oligomer poses against the saved Nyl12 PA6 register."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .chemistry import Substrate, atom_by_name


def pdb_coordinates(path: Path) -> dict[str, np.ndarray]:
    """Read named ligand coordinates from an existing PDB file.

    Args:
        path: Saved enzyme–substrate PDB containing a single LIG residue.

    Returns:
        Ligand atom names and their Cartesian coordinates.
    """
    coordinates: dict[str, np.ndarray] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[17:20] != "LIG":
            continue
        name = line[12:16].strip()
        if name in coordinates:
            raise ValueError(f"Duplicate ligand atom name {name} in {path}")
        coordinates[name] = np.array(
            [float(line[start : start + 8]) for start in (30, 38, 46)]
        )
    if not coordinates:
        raise ValueError(f"No LIG atoms found in {path}")
    return coordinates


def _protein_coordinates(path: Path) -> np.ndarray:
    """Read all protein heavy-atom coordinates from an existing PDB."""
    points = [
        [float(line[start : start + 8]) for start in (30, 38, 46)]
        for line in path.read_text().splitlines()
        if line.startswith("ATOM  ")
        and line[17:20] != "LIG"
        and line[76:78].strip().upper() not in {"H", "D"}
        and not line[12:16].strip().startswith("H")
    ]
    if not points:
        raise ValueError(f"No protein heavy atoms found in {path}")
    return np.asarray(points)


def _aligned_points(
    points: np.ndarray, source_anchor: np.ndarray, target_anchor: np.ndarray
) -> np.ndarray:
    """Align three ligand atoms to the bound PA6 amide by Kabsch rotation."""
    source_center = source_anchor.mean(axis=0)
    target_center = target_anchor.mean(axis=0)
    covariance = (source_anchor - source_center).T @ (target_anchor - target_center)
    u, _, vt = np.linalg.svd(covariance)
    handedness = np.linalg.det(u @ vt)
    rotation = u @ np.diag([1.0, 1.0, np.sign(handedness)]) @ vt
    return (points - source_center) @ rotation + target_center


def embed_on_pa6_register(
    substrate: Substrate,
    source_complex: Path,
    conformers: int = 24,
    seed: int = 233,
) -> tuple[Chem.Mol, dict[str, float | int]]:
    """Embed poly(X,AMA) and select the least clashing aligned conformer.

    The carbonyl C/O and leaving N are aligned to the existing PA6 L2
    producing register. This is an initial pose generator, not a binding score.

    Args:
        substrate: Named oligomer with explicit hydrogens.
        source_complex: Existing Nyl12–PA6 complex PDB.
        conformers: Number of independently embedded conformers.
        seed: Reproducible RDKit embedding seed.

    Returns:
        A copy with one aligned conformer and its clash audit.
    """
    if conformers < 1:
        raise ValueError("At least one conformer is required")
    mol = Chem.Mol(substrate.mol)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    parameters.numThreads = 4
    parameters.pruneRmsThresh = -1.0
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=conformers, params=parameters))
    if len(ids) < max(4, conformers // 2):
        raise ValueError(f"Only embedded {len(ids)} of {conformers} conformers")
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=4, maxIters=1000)
    if any(status < 0 for status, _ in results):
        raise ValueError("MMFF lacks parameters for one or more oligomer conformers")
    source_ligand = pdb_coordinates(source_complex)
    target_anchor = np.array([source_ligand[name] for name in ("C54", "O49", "N59")])
    source_indices = [
        atom_by_name(mol, name).GetIdx()
        for name in (substrate.carbonyl_c, substrate.carbonyl_o, substrate.leaving_n)
    ]
    protein = _protein_coordinates(source_complex)
    heavy_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    best: tuple[float, np.ndarray, int] | None = None
    for conformer_id in ids:
        points = np.asarray(mol.GetConformer(conformer_id).GetPositions())
        aligned = _aligned_points(points, points[source_indices], target_anchor)
        ligand_heavy = aligned[heavy_indices]
        distances = np.linalg.norm(ligand_heavy[:, None] - protein[None, :], axis=2)
        shortfall = np.maximum(2.4 - distances, 0.0)
        score = float(np.square(shortfall).sum())
        if best is None or score < best[0]:
            best = (score, aligned, conformer_id)
    if best is None:
        raise ValueError("No oligomer pose survived embedding")
    result = Chem.Mol(mol)
    result.RemoveAllConformers()
    conformer = Chem.Conformer(result.GetNumAtoms())
    for index, point in enumerate(best[1]):
        conformer.SetAtomPosition(index, point)
    result.AddConformer(conformer)
    return result, {
        "requested_conformers": conformers,
        "embedded_conformers": len(ids),
        "converged_conformers": sum(status == 0 for status, _ in results),
        "selected_mmff_converged": results[best[2]][0] == 0,
        "seed": seed,
        "selected_conformer": best[2],
        "protein_clash_penalty": best[0],
    }


def write_named_pdb(mol: Chem.Mol, path: Path) -> None:
    """Write named heavy-atom coordinates as one LIG residue.

    Args:
        mol: Molecule with one conformer and unique atom-name properties.
        path: Output PDB path.
    """
    conformer = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        name = atom.GetProp("name")
        point = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"HETATM{len(lines) + 1:5d} {name:4s} LIG A   1    "
            f"{point.x:8.3f}{point.y:8.3f}{point.z:8.3f}"
            f"  1.00  0.00          {atom.GetSymbol():>2s}"
        )
    path.write_text("\n".join([*lines, "END"]) + "\n")
