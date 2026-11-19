"""Generate matched DFTB3, GFN2-xTB, and DFT proton-transfer inputs."""

from __future__ import annotations

import argparse
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
DISTANCES = (2.5, 2.7, 3.0)
POINTS = 9


@dataclass(frozen=True)
class PairDefinition:
    """Chemical definition of one donor/acceptor benchmark pair."""

    name: str
    donor_smiles: str
    acceptor_smiles: str
    donor_anchor_atomic_number: int
    acceptor_anchor_atomic_number: int
    acceptor_formal_charge: int
    charge: int


PAIRS = (
    PairDefinition(
        "thr_to_terminal_amine",
        "CC(C)O",
        "CN",
        8,
        7,
        0,
        0,
    ),
    PairDefinition(
        "thr_to_asp",
        "CC(C)O",
        "CC(=O)[O-]",
        8,
        8,
        -1,
        -1,
    ),
    PairDefinition(
        "asph_to_terminal_amine",
        "CC(=O)O",
        "CN",
        8,
        7,
        0,
        0,
    ),
)


def _embed(smiles: str, seed: int) -> Chem.Mol:
    """Embed and UFF-optimize an explicit-hydrogen molecule."""
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise RuntimeError(f"RDKit embedding failed for {smiles}")
    if AllChem.UFFOptimizeMolecule(molecule, maxIters=1000) != 0:
        raise RuntimeError(f"UFF optimization failed for {smiles}")
    return molecule


def _anchor(molecule: Chem.Mol, atomic_number: int, formal_charge: int) -> int:
    """Return the unique atom matching an element and formal charge."""
    matches = [
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() == atomic_number
        and atom.GetFormalCharge() == formal_charge
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one Z={atomic_number}, charge={formal_charge}; got {matches}"
        )
    return matches[0]


def _donor_anchor(molecule: Chem.Mol, atomic_number: int) -> int:
    """Return the unique element-matched atom bonded to a hydrogen."""
    matches = [
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() == atomic_number
        and any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one protonated donor Z={atomic_number}; got {matches}"
        )
    return matches[0]


def _bonded_heavy(molecule: Chem.Mol, anchor: int) -> int:
    """Return the unique heavy atom directly bonded to an anchor atom."""
    matches = [
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(anchor).GetNeighbors()
        if neighbor.GetAtomicNum() > 1
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one heavy neighbor for atom {anchor}: {matches}")
    return matches[0]


def _bonded_hydrogen(molecule: Chem.Mol, anchor: int) -> int:
    """Return the unique hydrogen directly bonded to a donor anchor."""
    matches = [
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(anchor).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one donor hydrogen for atom {anchor}: {matches}")
    return matches[0]


def _rotation_between(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Return a proper rotation matrix aligning two nonzero vectors."""
    source = source / np.linalg.norm(source)
    target = target / np.linalg.norm(target)
    cross = np.cross(source, target)
    cosine = float(np.dot(source, target))
    if np.linalg.norm(cross) < 1.0e-10:
        if cosine > 0:
            return np.eye(3)
        axis = np.array([0.0, 1.0, 0.0])
        if abs(float(np.dot(source, axis))) > 0.9:
            axis = np.array([0.0, 0.0, 1.0])
        axis -= source * np.dot(source, axis)
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    skew = np.array(
        [[0.0, -cross[2], cross[1]], [cross[2], 0.0, -cross[0]], [-cross[1], cross[0], 0.0]]
    )
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / np.dot(cross, cross))


def _oriented_fragment(
    molecule: Chem.Mol,
    anchor: int,
    heavy_neighbor: int,
    target_anchor: np.ndarray,
    heavy_direction: np.ndarray,
) -> np.ndarray:
    """Translate and rotate a fragment around its reactive anchor."""
    coordinates = np.asarray(molecule.GetConformer().GetPositions(), dtype=float)
    centered = coordinates - coordinates[anchor]
    rotation = _rotation_between(centered[heavy_neighbor], heavy_direction)
    return centered @ rotation.T + target_anchor


def _symbols(molecule: Chem.Mol) -> list[str]:
    """Return atom symbols in RDKit atom order."""
    return [atom.GetSymbol() for atom in molecule.GetAtoms()]


def _render_xyz(symbols: list[str], coordinates: np.ndarray, title: str) -> str:
    """Render an XYZ coordinate file."""
    lines = [str(len(symbols)), title]
    lines.extend(
        f"{symbol:2s} {x: .10f} {y: .10f} {z: .10f}"
        for symbol, (x, y, z) in zip(symbols, coordinates, strict=True)
    )
    return "\n".join(lines) + "\n"


def _render_sqm(symbols: list[str], coordinates: np.ndarray, charge: int) -> str:
    """Render an Amber SQM DFTB3/3OB single-point input."""
    lines = [
        "Matched proton-transfer benchmark single point",
        " &qmmm",
        "  qm_theory='DFTB3',",
        f"  qmcharge={charge},",
        "  scfconv=1.0d-10, errconv=1.0d-6, tight_p_conv=1,",
        "  itrmax=500, maxcyc=0, verbosity=0,",
        " /",
    ]
    periodic = Chem.GetPeriodicTable()
    lines.extend(
        f" {periodic.GetAtomicNumber(symbol):2d} {symbol:2s} {x: .10f} {y: .10f} {z: .10f}"
        for symbol, (x, y, z) in zip(symbols, coordinates, strict=True)
    )
    return "\n".join(lines) + "\n"


def _render_gaussian(symbols: list[str], coordinates: np.ndarray, charge: int) -> str:
    """Render a Gaussian omegaB97X-D/aug-cc-pVTZ single-point input."""
    lines = [
        "%mem=2GB",
        "%nprocshared=2",
        "#p wB97XD/aug-cc-pVTZ SP SCF=(XQC,Tight) Integral=UltraFine NoSymm",
        "",
        "Matched proton-transfer benchmark single point",
        "",
        f"{charge} 1",
    ]
    lines.extend(
        f"{symbol:2s} {x: .10f} {y: .10f} {z: .10f}"
        for symbol, (x, y, z) in zip(symbols, coordinates, strict=True)
    )
    return "\n".join(lines) + "\n\n"


def generate(output: Path) -> None:
    """Generate all matched fixed-geometry benchmark inputs and metadata."""
    output.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, Any]] = []
    for pair_index, pair in enumerate(PAIRS, start=1):
        donor = _embed(pair.donor_smiles, 7100 + pair_index)
        acceptor = _embed(pair.acceptor_smiles, 7200 + pair_index)
        donor_anchor = _donor_anchor(donor, pair.donor_anchor_atomic_number)
        acceptor_anchor = _anchor(
            acceptor, pair.acceptor_anchor_atomic_number, pair.acceptor_formal_charge
        )
        donor_hydrogen = _bonded_hydrogen(donor, donor_anchor)
        donor_coordinates = _oriented_fragment(
            donor,
            donor_anchor,
            _bonded_heavy(donor, donor_anchor),
            np.zeros(3),
            np.array([-1.0, 0.0, 0.0]),
        )
        symbols = _symbols(donor) + _symbols(acceptor)
        for distance in DISTANCES:
            acceptor_coordinates = _oriented_fragment(
                acceptor,
                acceptor_anchor,
                _bonded_heavy(acceptor, acceptor_anchor),
                np.array([distance, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
            )
            for point, donor_h_distance in enumerate(
                np.linspace(0.95, distance - 1.02, POINTS)
            ):
                coordinates = np.vstack((donor_coordinates, acceptor_coordinates))
                coordinates[donor_hydrogen] = [donor_h_distance, 0.0, 0.0]
                identifier = f"{pair.name}_R{distance:.1f}_p{point:02d}"
                calculation = output / pair.name / f"R{distance:.1f}" / f"point_{point:02d}"
                calculation.mkdir(parents=True)
                title = f"{identifier}; O/N separation={distance:.3f} A"
                files = {
                    "xyz": calculation / "geometry.xyz",
                    "sqm": calculation / "dftb3.in",
                    "gaussian": calculation / "wb97xd_augccpvtz.gjf",
                    "charge": calculation / "charge.txt",
                }
                files["xyz"].write_text(_render_xyz(symbols, coordinates, title))
                files["sqm"].write_text(_render_sqm(symbols, coordinates, pair.charge))
                files["gaussian"].write_text(
                    _render_gaussian(symbols, coordinates, pair.charge)
                )
                files["charge"].write_text(f"{pair.charge}\n", encoding="utf-8")
                records.append(
                    {
                        "id": identifier,
                        "pair": pair.name,
                        "charge": pair.charge,
                        "donor_acceptor_A": distance,
                        "point": point,
                        "donor_H_A": float(donor_h_distance),
                        "acceptor_H_A": float(distance - donor_h_distance),
                        "directory": str(calculation.relative_to(output)),
                        "sha256": {
                            name: hashlib.sha256(path.read_bytes()).hexdigest()
                            for name, path in files.items()
                        },
                    }
                )
    metadata = {
        "purpose": "Fixed-geometry Hamiltonian comparison; not a free-energy profile",
        "systems": [pair.__dict__ for pair in PAIRS],
        "distances_A": DISTANCES,
        "points_per_distance": POINTS,
        "records": records,
    }
    (output / "benchmark_manifest.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info("Generated %d matched benchmark geometries in %s", len(records), output)


def main() -> None:
    """Parse the output path and generate the benchmark."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    generate(args.output)


if __name__ == "__main__":
    main()
