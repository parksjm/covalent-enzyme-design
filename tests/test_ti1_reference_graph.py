"""Tests for the coupled TI1 reference chemical graph."""

from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdchem

from nyl12_covalent.ti1_reference_graph import build_ti1_reference_graph

PA6_MOL2 = (
    Path(__file__).resolve().parents[1]
    / "provenance"
    / "xa_ma_capped_reference_20261022"
    / "pa6_l4_source.mol2"
)


def _origin(molecule: Chem.Mol, origin_id: str) -> Chem.Atom:
    matches = [
        atom for atom in molecule.GetAtoms() if atom.GetProp("origin_id") == origin_id
    ]
    assert len(matches) == 1
    return matches[0]


def test_builds_neutral_covalent_graph_with_stable_atom_origins() -> None:
    molecule, metadata = build_ti1_reference_graph(PA6_MOL2, "R")

    assert molecule.GetNumAtoms() == 103
    assert Chem.GetFormalCharge(molecule) == 0
    assert metadata["formula"] == "C29H60N6O8"
    assert metadata["attacked_carbon_cip"] == "R"
    assert len(metadata["atom_index_by_origin"]) == 103
    assert _origin(molecule, "THR:HG1").GetProp("atom_name") == "H3"
    assert _origin(molecule, "PA6:H3").GetIdx() != _origin(molecule, "THR:HG1").GetIdx()
    assert _origin(molecule, "PA6:H3").GetNeighbors()[0].GetProp("name") == "N56"


def test_connectivity_charges_stereochemistry_and_water_are_explicit() -> None:
    molecule, _ = build_ti1_reference_graph(PA6_MOL2, "S")
    atom = {entry.GetProp("origin_id"): entry for entry in molecule.GetAtoms()}
    assert atom["PA6:C54"].GetProp("_CIPCode") == "S"
    assert atom["THR:N"].GetFormalCharge() == 1
    assert atom["THR:OG1"].GetFormalCharge() == 0
    assert atom["PA6:O49"].GetFormalCharge() == -1
    assert (
        molecule.GetBondBetweenAtoms(
            atom["THR:OG1"].GetIdx(), atom["PA6:C54"].GetIdx()
        ).GetBondType()
        == rdchem.BondType.SINGLE
    )
    assert (
        molecule.GetBondBetweenAtoms(
            atom["PA6:C54"].GetIdx(), atom["PA6:O49"].GetIdx()
        ).GetBondType()
        == rdchem.BondType.SINGLE
    )
    assert (
        molecule.GetBondBetweenAtoms(
            atom["PA6:C54"].GetIdx(), atom["PA6:N59"].GetIdx()
        ).GetBondType()
        == rdchem.BondType.SINGLE
    )
    assert molecule.GetBondBetweenAtoms(
        atom["THR:N"].GetIdx(), atom["THR:HG1"].GetIdx()
    )
    assert (
        molecule.GetBondBetweenAtoms(atom["THR:OG1"].GetIdx(), atom["THR:HG1"].GetIdx())
        is None
    )
    assert len(atom["WAT:O"].GetNeighbors()) == 2
    assert all(neighbor.GetSymbol() == "H" for neighbor in atom["WAT:O"].GetNeighbors())
    assert molecule.GetNumConformers() == 0


@pytest.mark.parametrize("choice", ["r", "", None])
def test_requires_explicit_supported_c54_configuration(choice: object) -> None:
    with pytest.raises(ValueError, match="explicitly 'R' or 'S'"):
        build_ti1_reference_graph(PA6_MOL2, choice)  # type: ignore[arg-type]


def test_missing_pa6_source_fails_explicitly(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="PA6 MOL2 source does not exist"):
        build_ti1_reference_graph(tmp_path / "missing.mol2", "R")
