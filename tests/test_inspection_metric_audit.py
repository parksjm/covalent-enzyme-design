"""Separate precision loss from the actual reversible coordinate conversion."""

from copy import deepcopy
from typing import Any

import pytest

from nyl12_covalent.inspection_metric_audit import serialization_metric_audit


def evaluation() -> dict[str, Any]:
    rows = [
        {
            "site": site,
            "ligand": ligand,
            "carbonyl_c": f"C{amide}",
            "carbonyl_o": f"O{amide}",
            "leaving_n": f"N{amide}",
            "d1": 30.0,
            "d2": 30.0,
            "d3": 30.0,
            "a_bd": 110.0,
            "herschlag_phi": 104.25550754414714,
            "tier": "reject",
        }
        for site in ("AE", "BF", "CG", "DH")
        for ligand in "IJKL"
        for amide in range(3)
    ]
    return {"all_amide_site_rows": rows, "assigned": [rows[i] for i in (0, 15, 30, 45)]}


def test_large_quantization_diagnostic_not_conversion_error() -> None:
    raw = evaluation()
    pdb = deepcopy(raw)
    pdb["all_amide_site_rows"][0]["herschlag_phi"] = 103.60698284568883
    result = serialization_metric_audit(raw, pdb, deepcopy(pdb))
    assert result["metric_max_errors"]["herschlag_phi"] == 0
    assert result["native_to_PDB_quantization_metric_max_errors"][
        "herschlag_phi"
    ] == pytest.approx(0.6485246984582886)
    assert len(result["native_to_PDB_quantization_all48"]) == 48
    assert result["native_to_PDB_positive_AB_membership_changes"] == []


def test_actual_assigned_gate_changes_are_retained() -> None:
    raw = evaluation()
    raw["all_amide_site_rows"][0].update(tier="A", herschlag_phi=0.0001)
    pdb = deepcopy(raw)
    pdb["all_amide_site_rows"][0].update(tier="reject", herschlag_phi=-0.0001)
    pdb["assigned"] = [pdb["all_amide_site_rows"][i] for i in (1, 15, 30, 45)]
    result = serialization_metric_audit(raw, pdb, deepcopy(pdb))
    assert len(result["native_to_PDB_tier_changes"]) == 1
    changed = result["native_to_PDB_positive_AB_membership_changes"]
    assert len(changed) == 1 and changed[0]["native_legacy_assigned"]
    assert result["native_to_PDB_legacy_assignment_changed"]
    assert result["native_positive_phi_AB_gate"]["assigned_ab_count"] == 1
    assert result["serialized_positive_phi_AB_gate"]["assigned_ab_count"] == 0


@pytest.mark.parametrize(
    "metric,delta", [("d1", 0.0031), ("a_bd", 0.21), ("herschlag_phi", 0.21)]
)
def test_conversion_thresholds_unchanged(metric: str, delta: float) -> None:
    raw = evaluation()
    back = deepcopy(raw)
    back["all_amide_site_rows"][0][metric] += delta
    with pytest.raises(ValueError, match="Roundtrip metric mismatch"):
        serialization_metric_audit(raw, raw, back)


def test_circular_phi_and_linear_bd() -> None:
    raw = evaluation()
    raw["all_amide_site_rows"][0]["herschlag_phi"] = 179.99
    pdb = deepcopy(raw)
    pdb["all_amide_site_rows"][0]["herschlag_phi"] = -179.99
    result = serialization_metric_audit(raw, pdb, deepcopy(pdb))
    assert result["native_to_PDB_quantization_metric_max_errors"][
        "herschlag_phi"
    ] == pytest.approx(0.02)


def test_invalid_or_mismatched_rows_fail() -> None:
    raw = evaluation()
    bad = deepcopy(raw)
    bad["all_amide_site_rows"].pop()
    with pytest.raises(ValueError, match="48 unique"):
        serialization_metric_audit(raw, bad, raw)
    bad = deepcopy(raw)
    bad["all_amide_site_rows"][0]["carbonyl_c"] = "DIFFERENT"
    with pytest.raises(ValueError, match="identities differ"):
        serialization_metric_audit(raw, bad, raw)
