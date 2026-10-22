"""Tests for the local capped TI1 bonded-fitting surrogate."""

import pytest
from rdkit import Chem
from rdkit.Chem import rdchem

from nyl12_covalent.ti1_bonded_model import build_ti1_bonded_model


def _by_origin(molecule: Chem.Mol) -> dict[str, Chem.Atom]:
    return {atom.GetProp("origin_id"): atom for atom in molecule.GetAtoms()}


def test_model_has_expected_surrogate_inventory_and_scope() -> None:
    molecule, metadata = build_ti1_bonded_model("R")

    assert molecule.GetNumAtoms() == 37
    assert Chem.GetFormalCharge(molecule) == -1
    assert metadata["formula"] == "C10H20N3O4-"
    assert metadata["electron_count"] == 134
    assert metadata["multiplicity"] == 1
    assert metadata["model_scope"] == "bonded_fit_surrogate_only"
    assert metadata["native_esp_resp_transfer_eligible"] is False
    assert len(metadata["atom_index_by_origin"]) == 37
    assert molecule.GetNumConformers() == 0


def test_model_connectivity_caps_proton_reservoir_and_stereo_are_explicit() -> None:
    molecule, metadata = build_ti1_bonded_model("S")
    atom = _by_origin(molecule)
    bond = molecule.GetBondBetweenAtoms

    assert atom["THR:CA"].GetProp("_CIPCode") == "S"
    assert atom["THR:CB"].GetProp("_CIPCode") == "R"
    assert atom["TIL:C*"].GetProp("_CIPCode") == "S"
    assert atom["THR:N"].GetFormalCharge() == 0
    assert atom["THR:OG1"].GetFormalCharge() == 0
    assert atom["TIL:O"].GetFormalCharge() == -1
    assert "THR:HG1" not in atom
    assert metadata["proton_reservoir"]["omitted_origin_id"] == "THR:HG1"
    assert (
        bond(atom["THR:OG1"].GetIdx(), atom["TIL:C*"].GetIdx()).GetBondType()
        == rdchem.BondType.SINGLE
    )
    assert (
        bond(atom["TIL:C*"].GetIdx(), atom["TIL:O"].GetIdx()).GetBondType()
        == rdchem.BondType.SINGLE
    )
    assert (
        bond(atom["TIL:C*"].GetIdx(), atom["TIL:N"].GetIdx()).GetBondType()
        == rdchem.BondType.SINGLE
    )
    assert len(atom["THR:N"].GetNeighbors()) == 3
    assert len(atom["THR:OG1"].GetNeighbors()) == 2
    assert len(atom["TIL:N"].GetNeighbors()) == 3
    assert (
        bond(atom["ACE:C"].GetIdx(), atom["ACE:O"].GetIdx()).GetBondType()
        == rdchem.BondType.DOUBLE
    )
    assert (
        bond(atom["ACE:C"].GetIdx(), atom["THR:N"].GetIdx()).GetBondType()
        == rdchem.BondType.SINGLE
    )
    assert atom["NME:N"].GetProp("residue_id") == "NME"
    assert atom["THR:N"].GetProp("residue_id") == "THR"


def test_attacked_carbon_configuration_is_a_required_choice() -> None:
    r_model, r_metadata = build_ti1_bonded_model("R")
    s_model, s_metadata = build_ti1_bonded_model("S")

    assert _by_origin(r_model)["TIL:C*"].GetProp("_CIPCode") == "R"
    assert _by_origin(s_model)["TIL:C*"].GetProp("_CIPCode") == "S"
    assert r_metadata["attacked_carbon_cip"] == "R"
    assert s_metadata["attacked_carbon_cip"] == "S"


@pytest.mark.parametrize("choice", ["r", "", None])
def test_rejects_missing_or_unsupported_stereochemistry(choice: object) -> None:
    with pytest.raises(ValueError, match="explicitly 'R' or 'S'"):
        build_ti1_bonded_model(choice)  # type: ignore[arg-type]
