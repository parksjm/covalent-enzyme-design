"""Prepare explicit nylon-6 fragments for GAFF2 torsion benchmarking."""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

LOGGER = logging.getLogger(__name__)
RANDOM_SEED = 20261028


@dataclass(frozen=True)
class FragmentSpec:
    """Declare one chemically explicit nylon-6 benchmark fragment."""

    label: str
    smiles: str
    charge: int
    purpose: str
    proper_torsions: dict[str, tuple[int, int, int, int]]


FRAGMENTS = (
    FragmentSpec(
        label="internal_neutral",
        smiles="CC(=O)NCCCCCC(=O)NC",
        charge=0,
        purpose="Internal nylon-6 amides and methylene backbone",
        proper_torsions={
            "amide_omega": (0, 1, 3, 4),
            "amine_side": (1, 3, 4, 5),
            "central_methylene": (4, 5, 6, 7),
            "carbonyl_side": (6, 7, 8, 9),
        },
    ),
    FragmentSpec(
        label="n_terminus_cation",
        smiles="[NH3+]CCCCCC(=O)NC",
        charge=1,
        purpose="Solvent-exposed charged PA6 amino terminus",
        proper_torsions={
            "terminal_backbone": (0, 1, 2, 3),
            "ammonium_rotation": (2, 1, 0, 10),
        },
    ),
    FragmentSpec(
        label="c_terminus_anion",
        smiles="CC(=O)NCCCCCC(=O)[O-]",
        charge=-1,
        purpose="Solvent-exposed charged PA6 carboxylate terminus",
        proper_torsions={"terminal_backbone": (6, 7, 8, 9)},
    ),
)


def prepare_torsion_fragments(output_dir: Path) -> dict[str, Any]:
    """Generate deterministic MMFF seeds and Gaussian optimization inputs.

    Args:
        output_dir: New directory for fragment SDF/XYZ/Gaussian inputs.

    Returns:
        Serializable fragment manifest.

    Raises:
        FileExistsError: If the output directory already exists.
        RuntimeError: If embedding or MMFF94s minimization fails.
        ValueError: If a declared graph, charge, or torsion is invalid.
    """
    if output_dir.exists():
        raise FileExistsError(f"Fragment output already exists: {output_dir}")
    records = []
    output_dir.mkdir(parents=True)
    for spec in FRAGMENTS:
        molecule, energy = _embed_low_energy_conformer(spec)
        fragment_dir = output_dir / spec.label
        fragment_dir.mkdir()
        sdf_path = fragment_dir / "seed_mmff94s.sdf"
        xyz_path = fragment_dir / "seed_mmff94s.xyz"
        gaussian_path = fragment_dir / "opt_smd.gjf"
        with Chem.SDWriter(str(sdf_path)) as writer:
            writer.write(molecule)
        xyz_path.write_text(_xyz_text(molecule, spec.label))
        gaussian_path.write_text(_gaussian_input(molecule, spec))
        records.append(
            {
                "label": spec.label,
                "smiles": spec.smiles,
                "charge": spec.charge,
                "multiplicity": 1,
                "purpose": spec.purpose,
                "atom_count": molecule.GetNumAtoms(),
                "atom_names_in_order": [
                    atom.GetProp("atom_name") for atom in molecule.GetAtoms()
                ],
                "proper_torsions_zero_based": {
                    name: list(indices)
                    for name, indices in spec.proper_torsions.items()
                },
                "proper_torsions_gaussian_one_based": {
                    name: [index + 1 for index in indices]
                    for name, indices in spec.proper_torsions.items()
                },
                "seed_mmff94s_energy": energy,
                "sha256": {
                    "sdf": _sha256(sdf_path),
                    "xyz": _sha256(xyz_path),
                    "gaussian_input": _sha256(gaussian_path),
                },
            }
        )
    manifest = {
        "schema_version": 1,
        "scope": "PA6 L4 GAFF2 torsion/improper QM benchmark fragments",
        "random_seed": RANDOM_SEED,
        "seed_method": "ETKDGv3, 30 conformers, converged MMFF94s minimum",
        "qm_geometry_method": "wB97XD/def2SVP/SMD(water)",
        "fragments": records,
        "scientific_limits": [
            "Fragment scans validate local bonded terms, not full-oligomer folding.",
            "RESP charges must be generated separately for MM profile comparison.",
            "Bidirectional relaxed scans and continuity checks are required.",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared %d PA6 torsion fragments", len(records))
    return manifest


def _embed_low_energy_conformer(spec: FragmentSpec) -> tuple[Chem.Mol, float]:
    molecule = Chem.MolFromSmiles(spec.smiles)
    if molecule is None:
        raise ValueError(f"Invalid fragment SMILES: {spec.smiles}")
    molecule = Chem.AddHs(molecule)
    if Chem.GetFormalCharge(molecule) != spec.charge:
        raise ValueError(f"Unexpected formal charge for {spec.label}")
    for name, indices in spec.proper_torsions.items():
        if len(set(indices)) != 4 or any(
            molecule.GetBondBetweenAtoms(first, second) is None
            for first, second in zip(indices, indices[1:])
        ):
            raise ValueError(f"Invalid proper torsion {spec.label}:{name}")
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    conformers = list(AllChem.EmbedMultipleConfs(molecule, numConfs=30, params=params))
    if len(conformers) != 30 or not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise RuntimeError(f"Could not parameterize 30 seeds for {spec.label}")
    results = AllChem.MMFFOptimizeMoleculeConfs(
        molecule, mmffVariant="MMFF94s", maxIters=5000, numThreads=1
    )
    converged = [
        (energy, conformer_id)
        for conformer_id, (status, energy) in zip(conformers, results, strict=True)
        if status == 0 and np.isfinite(energy)
    ]
    if not converged:
        raise RuntimeError(f"No converged MMFF94s seed for {spec.label}")
    energy, conformer_id = min(converged)
    selected = Chem.Mol(molecule)
    selected.RemoveAllConformers()
    selected.AddConformer(molecule.GetConformer(conformer_id), assignId=True)
    _name_atoms(selected)
    return selected, float(energy)


def _name_atoms(molecule: Chem.Mol) -> None:
    counts: dict[str, int] = {}
    for atom in molecule.GetAtoms():
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
        atom.SetProp("atom_name", f"{symbol}{counts[symbol]:02d}")


def _xyz_text(molecule: Chem.Mol, label: str) -> str:
    coordinates = molecule.GetConformer().GetPositions()
    rows = [
        f"{atom.GetSymbol():2s} {x: .10f} {y: .10f} {z: .10f}"
        for atom, (x, y, z) in zip(molecule.GetAtoms(), coordinates, strict=True)
    ]
    return "\n".join([str(molecule.GetNumAtoms()), label, *rows, ""])


def _gaussian_input(molecule: Chem.Mol, spec: FragmentSpec) -> str:
    coordinates = molecule.GetConformer().GetPositions()
    lines = [
        "%nprocshared=8",
        "%mem=16GB",
        f"%chk={spec.label}_opt_smd.chk",
        "#p wB97XD/def2SVP SCRF=(SMD,Solvent=Water) "
        "Opt=(Tight,MaxCycles=300) SCF=(Tight,XQC) NoSymm",
        "",
        f"PA6 torsion benchmark fragment: {spec.label}",
        "",
        f"{spec.charge} 1",
    ]
    lines.extend(
        f"{atom.GetSymbol():2s} {x: .10f} {y: .10f} {z: .10f}"
        for atom, (x, y, z) in zip(molecule.GetAtoms(), coordinates, strict=True)
    )
    return "\n".join([*lines, "", ""])


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
