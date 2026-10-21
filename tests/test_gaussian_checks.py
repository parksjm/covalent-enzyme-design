"""Synthetic regression checks for Gaussian scientific acceptance gates."""

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from nyl12_covalent.gaussian_checks import (
    assess_gaussian,
    final_orientation,
    validate_geometry,
)


def _water(tmp_path: Path) -> tuple[Path, Chem.Mol, np.ndarray]:
    molecule = Chem.AddHs(Chem.MolFromSmiles("O"))
    xyz = np.array([[0.0, 0.0, 0.0], [0.957, 0.0, 0.0], [-0.24, 0.927, 0.0]])
    conf = Chem.Conformer(3)
    for index, row in enumerate(xyz):
        conf.SetAtomPosition(index, tuple(row))
    molecule.AddConformer(conf)
    path = tmp_path / "molecule.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(molecule)
    return path, molecule, xyz


def _log(xyz: np.ndarray) -> str:
    rows = "\n".join(
        f"{i + 1} {z} 0 {x} {y} {v}"
        for i, (z, (x, y, v)) in enumerate(zip([8, 1, 1], xyz))
    )
    return f"""Charge = 0 Multiplicity = 1
Input orientation:
 ------------------
 Center Atomic Type X Y Z
 ------------------
{rows}
 ------------------
Optimization completed.
Frequencies -- 1500.0 3500.0 3600.0
Normal termination of Gaussian 16
Electrostatic Properties (Atomic Units)
ESP Fit Center 1 is at 0.0 0.0 0.0
Normal termination of Gaussian 16
"""


def test_accepts_completed_minimum_and_preserved_graph(tmp_path: Path) -> None:
    sdf, _, xyz = _water(tmp_path)
    path = tmp_path / "g16.log"
    path.write_text(_log(xyz))
    report = assess_gaussian(path, sdf)
    assert report["accepted"]
    assert report["geometry"]["accepted"]
    assert report["final_coordinates_A"] == xyz.tolist()


@pytest.mark.parametrize(
    ("old", "new"),
    [
        ("1500.0", "-10.0"),
        ("1500.0 3500.0 3600.0", "3500.0 3600.0"),
        ("Optimization completed.", "Optimization stopped."),
        ("Multiplicity = 1", "Multiplicity = 3"),
        ("Charge = 0", "Charge = 1"),
        ("ESP Fit Center", "No ESP points"),
        ("Normal termination of Gaussian 16", "Error termination"),
    ],
)
def test_rejects_process_success_without_scientific_acceptance(
    tmp_path: Path, old: str, new: str
) -> None:
    sdf, _, xyz = _water(tmp_path)
    path = tmp_path / "g16.log"
    path.write_text(_log(xyz).replace(old, new))
    assert not assess_gaussian(path, sdf)["accepted"]


def test_checks_final_orientation_not_earlier_geometry(tmp_path: Path) -> None:
    _, mol, xyz = _water(tmp_path)
    moved = xyz.copy()
    moved[1, 0] = 4
    numbers, final = final_orientation(_log(xyz) + _log(moved))
    assert not validate_geometry(mol, numbers, final)["accepted"]


def test_rejects_new_short_contact(tmp_path: Path) -> None:
    _, mol, xyz = _water(tmp_path)
    xyz[2] = xyz[1] + [0.1, 0, 0]
    assert not validate_geometry(mol, [8, 1, 1], xyz)["accepted"]


@pytest.mark.parametrize(
    "text",
    [
        "missing",
        "Input orientation:\n",
        "Input orientation:\n-----\nheader\n-----\n2 8 0 0 0 0\n-----",
        "Input orientation:\n-----\nheader\n-----\n1 8 0 nan 0 0\n-----",
    ],
)
def test_rejects_missing_or_bad_coordinate_tables(text: str) -> None:
    with pytest.raises(ValueError):
        final_orientation(text)


def test_rejects_changed_inventory_and_nonfinite_coordinates(tmp_path: Path) -> None:
    _, mol, xyz = _water(tmp_path)
    with pytest.raises(ValueError, match="inventory"):
        validate_geometry(mol, [8, 1, 6], xyz)
    with pytest.raises(ValueError, match="shape"):
        validate_geometry(mol, [8, 1, 1], np.full((3, 3), np.nan))


def test_rejects_malformed_sdf(tmp_path: Path) -> None:
    sdf = tmp_path / "bad.sdf"
    sdf.write_text("not a molecule\n$$$$\n")
    with pytest.raises(ValueError, match="one valid"):
        assess_gaussian(tmp_path / "unused.log", sdf)


def test_detects_stereoinversion() -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)F"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    xyz = mol.GetConformer().GetPositions()
    numbers = [a.GetAtomicNum() for a in mol.GetAtoms()]
    assert validate_geometry(mol, numbers, xyz)["accepted"]
    xyz[:, 0] *= -1
    assert not validate_geometry(mol, numbers, xyz)["accepted"]
