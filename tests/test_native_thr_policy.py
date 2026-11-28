"""Regress the distinction between catalytic eligibility and native diagnostics."""

from __future__ import annotations

import pytest

from nyl12_covalent.native_thr_policy import geometry_eligible, policy_metadata


@pytest.mark.parametrize("tier", ["A", "B"])
def test_catalytic_geometry_is_eligible(tier: str) -> None:
    """Eligibility depends on catalytic geometry, not native similarity."""
    assert geometry_eligible(tier)


@pytest.mark.parametrize("tier", ["C", "none", "reject"])
def test_geometry_rejections_still_apply(tier: str) -> None:
    """Diagnostic-only native treatment does not accept bad catalytic geometry."""
    assert not geometry_eligible(tier)


@pytest.mark.parametrize("tier", ["", "unknown", None])
def test_invalid_tier_fails(tier: str) -> None:
    """Missing scientific classification cannot silently become eligible."""
    with pytest.raises(ValueError, match="Unrecognized"):
        geometry_eligible(tier)


def test_policy_is_explicit_and_metadata_cannot_mutate_policy() -> None:
    """Native thresholds remain diagnostic and unchanged in independent reports."""
    first = policy_metadata()
    assert first["native_thr_role"] == "diagnostic_only"
    assert first["native_thresholds_changed"] is False
    first["geometry_eligible_tiers"].clear()
    assert policy_metadata()["geometry_eligible_tiers"] == ["A", "B"]
