"""Calculate protein residue SASA and gate design-mask candidates."""

from __future__ import annotations

import logging
from pathlib import Path
import re
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFreeSASA

LOGGER = logging.getLogger(__name__)

TIEN2013_MAX_ASA: dict[str, float] = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLU": 223.0,
    "GLN": 225.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}


def protein_molecule_from_pdb(path: Path) -> Chem.Mol:
    """Load protein ATOM records with coordinates and atom metadata only.

    Args:
        path: PDB path containing protein atoms and possibly HETATM records.

    Returns:
        RDKit molecule with conformer and zero inferred bonds.

    Raises:
        ValueError: If no protein atoms or nonstandard residues are present.
    """
    lines = [
        line
        for line in path.read_text().splitlines()
        if line.startswith("ATOM  ") and line[76:78].strip().upper() not in {"H", "D"}
    ]
    if not lines:
        raise ValueError(f"No protein ATOM records in {path}")
    block = "\n".join(lines) + "\nEND\n"
    mol = Chem.MolFromPDBBlock(
        block, sanitize=False, removeHs=False, proximityBonding=False
    )
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"RDKit could not parse protein atoms from {path}")
    if mol.GetNumBonds() != 0:
        raise ValueError("PDB parsing unexpectedly created bonds")
    residues: set[tuple[str, int, str]] = set()
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError(f"Missing PDB atom metadata in {path}")
        residue = info.GetResidueName().strip()
        if residue not in TIEN2013_MAX_ASA:
            raise ValueError(f"Nonstandard amino acid {residue!r} in {path}")
        residues.add((info.GetChainId().strip(), info.GetResidueNumber(), residue))
    LOGGER.info("Loaded %d protein residues from %s", len(residues), path)
    return mol


def calculate_residue_sasa(mol: Chem.Mol) -> dict[tuple[str, int, str], float]:
    """Calculate residue ASA using Lee–Richards and Protor-classified radii.

    Args:
        mol: Coordinate-bearing protein molecule parsed without inferred bonds.

    Returns:
        Total solvent accessible surface area per chain/residue identity.

    Raises:
        ValueError: If classification fails or any radius/area is invalid.
    """
    options = rdFreeSASA.SASAOpts()
    options.algorithm = rdFreeSASA.SASAAlgorithm.LeeRichards
    options.classifier = rdFreeSASA.SASAClassifier.Protor
    options.probeRadius = 1.4
    radii = rdFreeSASA.classifyAtoms(mol, options)
    if radii is None or len(radii) != mol.GetNumAtoms():
        raise ValueError("Protor radii classification failed")
    if any(radius is None or radius <= 0 for radius in radii):
        raise ValueError("Protor returned a missing or nonpositive atomic radius")
    rdFreeSASA.CalcSASA(mol, radii, opts=options)
    residue_asa: dict[tuple[str, int, str], float] = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        assert info is not None
        key = (
            info.GetChainId().strip(),
            info.GetResidueNumber(),
            info.GetResidueName().strip(),
        )
        atom_asa = atom.GetDoubleProp("SASA")
        if atom_asa < 0:
            raise ValueError(f"Negative atomic SASA for {key}")
        residue_asa[key] = residue_asa.get(key, 0.0) + atom_asa
    return residue_asa


def relative_accessibility(residue: str, asa: float) -> float:
    """Normalize residue ASA by the Tien et al. allowed maximum ASA table.

    Args:
        residue: Three-letter standard amino-acid code.
        asa: Residue total solvent accessible surface area in square angstroms.

    Returns:
        Relative solvent accessibility.

    Raises:
        ValueError: If residue is unsupported or ASA is negative.
    """
    if residue not in TIEN2013_MAX_ASA or asa < 0:
        raise ValueError(f"Invalid residue/ASA pair: {residue!r}, {asa}")
    return asa / TIEN2013_MAX_ASA[residue]


def passes_all_copy_sasa(copy_rsas: list[float], threshold: float = 0.25) -> bool:
    """Return whether every equivalent-copy RSA meets the conservative gate.

    Args:
        copy_rsas: Relative accessibilities for all copies and parent structures.
        threshold: Minimum required RSA for every copy.

    Returns:
        True only when all supplied copies meet the threshold.

    Raises:
        ValueError: If the evidence list is empty or contains invalid values.
    """
    if not copy_rsas or not 0 <= threshold <= 1:
        raise ValueError("RSA evidence must be nonempty and threshold in [0, 1]")
    if any(not 0 <= value <= 2 for value in copy_rsas):
        raise ValueError("RSA evidence contains an invalid value")
    return all(value >= threshold for value in copy_rsas)


def read_bound_complex_heavy_coordinates(
    path: Path,
) -> tuple[dict[str, dict[int, tuple[str, np.ndarray]]], dict[str, np.ndarray]]:
    """Read protein and bound-substrate heavy coordinates from an eight-copy PDB.

    Args:
        path: PDB containing protein chains A–H and ligand chains I–L.

    Returns:
        Protein coordinates keyed by chain and residue, and ligand coordinates
        keyed by chain.

    Raises:
        ValueError: If expected protein or ligand chains are missing.
    """
    protein: dict[str, dict[int, list[tuple[str, np.ndarray]]]] = {}
    ligand: dict[str, list[np.ndarray]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[16] not in (" ", "A"):
            continue
        element = (
            line[76:78].strip().upper()
            or line[12:16].strip().upper().lstrip("0123456789")[:1]
        )
        if element in {"H", "D"}:
            continue
        chain, number, residue = line[21], int(line[22:26]), line[17:20].strip()
        xyz = np.array([float(line[i : i + 8]) for i in (30, 38, 46)])
        if line.startswith("ATOM  ") and chain in "ABCDEFGH":
            protein.setdefault(chain, {}).setdefault(number, []).append((residue, xyz))
        elif line.startswith("HETATM") and chain in "IJKL":
            ligand.setdefault(chain, []).append(xyz)
    if set(protein) != set("ABCDEFGH") or set(ligand) != set("IJKL"):
        raise ValueError(f"Missing expected chains in {path}")
    normalized = {
        chain: {
            number: (atoms[0][0], np.array([xyz for _, xyz in atoms]))
            for number, atoms in residues.items()
        }
        for chain, residues in protein.items()
    }
    return normalized, {chain: np.array(coords) for chain, coords in ligand.items()}


def proposal_shell_positions(
    proposal: dict[str, Any], seed: int, family: str, shell: str
) -> set[int]:
    """Read candidate residue positions from the parent/mask proposal JSON.

    Args:
        proposal: Parsed parent and mask proposal document.
        seed: WT model seed identifying the shell record.
        family: Protein family key, either ``alpha`` or ``beta``.
        shell: Distance shell, ``6-12`` or ``12-16``.

    Returns:
        Residue numbers proposed in the requested shell.

    Raises:
        ValueError: If the seed or shell is unsupported.
    """
    seed_key = {2026112202: "seed_02", 2026112203: "seed_03"}.get(seed)
    if seed_key is None or shell not in {"6-12", "12-16"}:
        raise ValueError(f"Unsupported seed/shell: {seed}, {shell}")
    record = proposal["mask_proposal"]["shell_geometry_and_candidates"][family][
        seed_key
    ]
    labels = (
        [row["residue"] for row in record["candidate_6_12_positions"]]
        if shell == "6-12"
        else record["distant_12_16_unselected"]
    )
    return {int(re.search(r"\d+$", label).group()) for label in labels}


def summarize_sasa(
    mol: Chem.Mol, expected_residues: int | None = None
) -> dict[str, Any]:
    """Return auditable residue ASA and RSA records for a protein molecule.

    Args:
        mol: Protein molecule with PDB atom metadata and calculated coordinates.
        expected_residues: Optional required number of distinct residues.

    Returns:
        JSON-compatible settings and per-residue ASA/RSA records.

    Raises:
        ValueError: If the distinct residue count differs from expectation.
    """
    residue_asa = calculate_residue_sasa(mol)
    if expected_residues is not None and len(residue_asa) != expected_residues:
        raise ValueError(
            f"Expected {expected_residues} standard residues, found {len(residue_asa)}"
        )
    records = [
        {
            "chain": chain,
            "residue_number": number,
            "residue_name": residue,
            "asa_A2": round(asa, 4),
            "rsa": round(relative_accessibility(residue, asa), 6),
        }
        for (chain, number, residue), asa in sorted(residue_asa.items())
    ]
    return {
        "algorithm": "Lee-Richards",
        "classifier": "Protor",
        "probe_radius_A": 1.4,
        "protein_records": "PDB ATOM records only; HETATM ligands excluded",
        "hydrogen_records": "PDB H and D records excluded",
        "reference": "Tien et al. 2013, PLOS ONE 8:e80635, Table 1",
        "residue_count": len(records),
        "residues": records,
    }
