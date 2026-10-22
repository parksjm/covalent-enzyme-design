"""Checks for the TI1 surrogate/native atom crosswalk."""

from pathlib import Path

import pytest
from rdkit import Chem

from nyl12_covalent.ti1_bonded_model import build_ti1_bonded_model
from nyl12_covalent.ti1_model_crosswalk import build_ti1_model_crosswalk
from nyl12_covalent.ti1_reference_graph import build_ti1_reference_graph


SOURCE = (
    Path(__file__).parents[1]
    / "provenance"
    / "xa_ma_capped_reference_20261022"
    / "pa6_l4_source.mol2"
)


def _graphs(stereochemistry: str) -> tuple[Chem.Mol, Chem.Mol]:
    """Build the bonded surrogate and source-derived native graph."""
    surrogate, _ = build_ti1_bonded_model(stereochemistry)
    native, _ = build_ti1_reference_graph(SOURCE, stereochemistry)
    return surrogate, native


@pytest.mark.parametrize("stereochemistry", ["R", "S"])
def test_crosswalk_reports_intended_correspondences(stereochemistry: str) -> None:
    """Both stereoisomers map heavy atoms while keeping charges ineligible."""
    surrogate, native = _graphs(stereochemistry)
    report = build_ti1_model_crosswalk(surrogate, native)

    assert len(report["correspondences"]) == 14
    assert report["native_charge_transfer_eligible"] is False
    assert report["forcefield_validated"] is False
    assert report["native_coordinates_generated"] is False
    assert report["surrogate_graph"] == {
        "atom_count": 37,
        "heavy_atom_count": 17,
        "formal_charge": -1,
    }
    assert report["native_graph"] == {
        "atom_count": 103,
        "heavy_atom_count": 43,
        "formal_charge": 0,
    }
    rows = {row["surrogate_origin_id"]: row for row in report["correspondences"]}
    assert rows["TIL:C*"]["native_origin_id"] == "PA6:C54"
    assert rows["TIL:leavingCH3"]["native_origin_id"] != "PA6:C54"
    assert all(row["native_charge_transfer_eligible"] is False for row in rows.values())
    assert all(row["forcefield_validated"] is False for row in rows.values())
    assert "THR:N protonation" in report["state_differences"][0]


def test_duplicate_origin_id_is_rejected() -> None:
    """The source identity contract rejects ambiguous atom origins."""
    surrogate, native = _graphs("R")
    duplicate = Chem.Mol(surrogate)
    duplicate.GetAtomWithIdx(1).SetProp(
        "origin_id", duplicate.GetAtomWithIdx(0).GetProp("origin_id")
    )

    with pytest.raises(ValueError, match="Duplicate surrogate origin_id"):
        build_ti1_model_crosswalk(duplicate, native)


def test_wrong_covalent_bond_is_rejected() -> None:
    """A broken reactive center cannot produce an apparently valid crosswalk."""
    surrogate, native = _graphs("R")
    broken = Chem.RWMol(surrogate)
    origin_indices = {
        atom.GetProp("origin_id"): atom.GetIdx() for atom in surrogate.GetAtoms()
    }
    broken.RemoveBond(origin_indices["TIL:C*"], origin_indices["TIL:O"])

    with pytest.raises(ValueError, match=r"Surrogate bond TIL:C\*-TIL:O"):
        build_ti1_model_crosswalk(broken.GetMol(), native)
