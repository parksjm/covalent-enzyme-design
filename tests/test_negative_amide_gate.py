"""Behavioral checks for chemically mapped negative-phi pose assignment."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import pytest

from nyl12_covalent.negative_amide_gate import evaluate_negative_amide_gate

SITES = ("AE", "BF")
LIGANDS = ("I", "J")


def make_rows(
    sites: tuple[str, ...] = SITES, ligands: tuple[str, ...] = LIGANDS
) -> list[dict[str, Any]]:
    """Create every physical site/ligand/internal-amide combination."""
    return [
        {
            "site": site,
            "ligand": ligand,
            "carbonyl_c": f"C{amide}",
            "carbonyl_o": f"O{amide}",
            "leaving_n": f"N{amide}",
            "herschlag_phi": 85.0,
            "d1": 3.0,
            "tier": "reject",
        }
        for site in sites
        for ligand in ligands
        for amide in range(3)
    ]


def set_pose(
    rows: list[dict[str, Any]],
    site: str,
    ligand: str,
    amide: int = 0,
    *,
    tier: str = "A",
    phi: float = -85.0,
    d1: float = 3.0,
) -> None:
    """Set one explicitly identified pose in a complete geometry fixture."""
    row = next(
        row
        for row in rows
        if row["site"] == site
        and row["ligand"] == ligand
        and row["carbonyl_c"] == f"C{amide}"
    )
    row.update(tier=tier, herschlag_phi=phi, d1=d1)


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Evaluate the compact complete two-by-two physical system."""
    return evaluate_negative_amide_gate(rows, sites=SITES, ligands=LIGANDS)


def test_default_48_rows_reject_positive_only_legacy_a() -> None:
    rows = make_rows(("AE", "BF", "CG", "DH"), ("I", "J", "K", "L"))
    for row in rows:
        row["tier"] = "A"
    result = evaluate_negative_amide_gate(rows)
    assert result["input_row_count"] == 48
    assert result["eligible_negative_ab_count"] == 0
    assert result["assigned_negative_ab_count"] == 0
    assert result["assigned"] == [None] * 4
    assert result["passes_pose_gate"] is False


def test_filter_before_pair_ranking_retains_hidden_negative_b() -> None:
    rows = make_rows()
    set_pose(rows, "AE", "I", tier="A", phi=85.0, d1=2.5)
    set_pose(rows, "AE", "I", 1, tier="B", d1=3.1)
    set_pose(rows, "AE", "I", 2, tier="C", d1=2.4)
    set_pose(rows, "BF", "J", tier="A", phi=0.0)
    result = evaluate(rows)
    assert result["eligible_negative_ab_count"] == 1
    assert result["assigned"][0]["carbonyl_c"] == "C1"
    assert result["assigned"][0]["tier"] == "B"
    assert result["assigned"][1] is None
    assert result["passes_pose_gate"] is True


def test_global_cardinality_beats_shortest_distance_greedy() -> None:
    rows = make_rows()
    set_pose(rows, "AE", "I", d1=2.5)
    set_pose(rows, "AE", "J", tier="B", d1=3.2)
    set_pose(rows, "BF", "I", tier="B", d1=3.3)
    result = evaluate(rows)
    assert result["assigned_negative_ab_count"] == 2
    assert [row["ligand"] for row in result["assigned"]] == ["J", "I"]
    assert result["assigned_site_order"] == list(SITES)


def test_tier_priority_then_total_distance_select_global_assignment() -> None:
    rows = make_rows()
    set_pose(rows, "AE", "I", d1=4.0)
    set_pose(rows, "BF", "J", d1=4.0)
    set_pose(rows, "AE", "J", tier="B", d1=2.0)
    set_pose(rows, "BF", "I", d1=2.0)
    assert [r["ligand"] for r in evaluate(rows)["assigned"]] == ["I", "J"]
    set_pose(rows, "AE", "J", tier="A", d1=2.0)
    assert [r["ligand"] for r in evaluate(rows)["assigned"]] == ["J", "I"]


def test_pair_ranking_tier_distance_then_atom_identity() -> None:
    rows = make_rows()
    set_pose(rows, "AE", "I", 0, tier="B", d1=2.5)
    set_pose(rows, "AE", "I", 1, tier="A", d1=3.2)
    set_pose(rows, "AE", "I", 2, tier="A", d1=3.0)
    assert evaluate(rows)["assigned"][0]["carbonyl_c"] == "C2"
    set_pose(rows, "AE", "I", 1, tier="A", d1=3.0)
    assert evaluate(rows)["assigned"][0]["carbonyl_c"] == "C1"


def test_source_unmodified_and_assignment_deterministic_under_reordering() -> None:
    rows = make_rows()
    for row in rows:
        row.update(tier="A", herschlag_phi=-85.0)
    original = deepcopy(rows)
    result = evaluate(rows)
    assert result == evaluate(list(reversed(rows)))
    assert [row["ligand"] for row in result["assigned"]] == ["I", "J"]
    assert rows == original
    assert result["legacy_tiers_modified"] is False
    result["assigned"][0]["tier"] = "reject"
    assert rows == original


@pytest.mark.parametrize(
    "field,value",
    [
        ("tier", "unknown"),
        ("herschlag_phi", float("nan")),
        ("herschlag_phi", float("inf")),
        ("d1", float("nan")),
        ("d1", float("inf")),
        ("d1", 0.0),
        ("d1", -1.0),
        ("site", "XX"),
        ("ligand", "X"),
        ("carbonyl_c", ""),
    ],
)
def test_invalid_row_fails_explicitly(field: str, value: Any) -> None:
    rows = make_rows()
    rows[0][field] = value
    with pytest.raises(ValueError):
        evaluate(rows)


def test_incomplete_rows_fail() -> None:
    with pytest.raises(ValueError, match="Expected all"):
        evaluate(make_rows()[:-1])


def test_duplicate_amide_rows_fail() -> None:
    rows = make_rows()
    rows[1] = dict(rows[0])
    with pytest.raises(ValueError, match="duplicate"):
        evaluate(rows)


def test_inconsistent_graph_amide_across_sites_fails() -> None:
    rows = make_rows()
    rows[0]["carbonyl_c"] = "C999"
    with pytest.raises(ValueError, match="Inconsistent graph-amide"):
        evaluate(rows)


def test_incomplete_pair_coverage_with_correct_total_fails() -> None:
    rows = make_rows()
    rows[0].update(site="BF", carbonyl_c="C999")
    with pytest.raises(ValueError, match="Incomplete site/ligand coverage"):
        evaluate(rows)


@pytest.mark.parametrize(
    "sites,ligands,amides",
    [
        ((), (), 3),
        (("AE", "AE"), LIGANDS, 3),
        (SITES, ("I", "I"), 3),
        (SITES, ("I",), 3),
        (SITES, LIGANDS, 0),
    ],
)
def test_invalid_physical_inventory_fails(
    sites: tuple[str, ...], ligands: tuple[str, ...], amides: int
) -> None:
    with pytest.raises(ValueError, match="unique, equally sized"):
        evaluate_negative_amide_gate(
            [], sites=sites, ligands=ligands, amides_per_ligand=amides
        )


@pytest.mark.parametrize("sign", [0, 2, -2, True, False, 1.0, -1.0, "1", None])
def test_signed_gate_requires_exact_explicit_integer_sign(sign: Any) -> None:
    from nyl12_covalent.negative_amide_gate import evaluate_signed_amide_gate

    with pytest.raises(ValueError, match="target_phi_sign"):
        evaluate_signed_amide_gate(
            make_rows(), target_phi_sign=sign, sites=SITES, ligands=LIGANDS
        )


@pytest.mark.parametrize("phi,expected", [(85.0, True), (-85.0, False), (0.0, False)])
def test_positive_campaign_accepts_only_actual_positive_phi(
    phi: float, expected: bool
) -> None:
    from nyl12_covalent.negative_amide_gate import evaluate_signed_amide_gate

    rows = make_rows()
    set_pose(rows, "AE", "I", phi=phi)
    original = deepcopy(rows)
    result = evaluate_signed_amide_gate(
        rows, target_phi_sign=+1, sites=SITES, ligands=LIGANDS
    )
    assert result["target_phi_sign"] == +1
    assert result["passes_pose_gate"] is expected
    assert result["assigned_ab_count"] == int(expected)
    assert rows == original
    if expected:
        assert result["eligible_ab"][0]["herschlag_phi"] == 85.0
        result["assigned"][0]["herschlag_phi"] = -999.0
        assert rows == original


def test_positive_hidden_b_survives_opposite_sign_a_ranking() -> None:
    from nyl12_covalent.negative_amide_gate import evaluate_signed_amide_gate

    rows = make_rows()
    set_pose(rows, "AE", "I", tier="A", phi=-88.0, d1=2.5)
    set_pose(rows, "AE", "I", 1, tier="B", phi=88.228, d1=3.1)
    result = evaluate_signed_amide_gate(
        rows, target_phi_sign=+1, sites=SITES, ligands=LIGANDS
    )
    assert result["eligible_ab_count"] == 1
    assert result["assigned_ab_count"] == 1
    assert result["assigned"][0]["carbonyl_c"] == "C1"
    assert result["assigned"][0]["herschlag_phi"] == 88.228
    assert result["assigned"][0]["tier"] == "B"
