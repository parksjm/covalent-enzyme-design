"""Separate catalytic-geometry eligibility from crystal-similarity diagnostics."""

from __future__ import annotations

from typing import Any

NATIVE_THR_ROLE = "diagnostic_only"


def geometry_eligible(tier: str) -> bool:
    """Return original Tier A/B eligibility without a native-Thr veto.

    Args:
        tier: Explicit original geometry classification, including reject labels.

    Returns:
        Whether the original catalytic geometry belongs to Tier A or B.

    Raises:
        ValueError: The supplied classification is missing or unrecognized.
    """
    if tier not in {"A", "B", "C", "reject", "none"}:
        raise ValueError(f"Unrecognized catalytic geometry tier: {tier!r}")
    return tier in {"A", "B"}


def policy_metadata() -> dict[str, Any]:
    """Return explicit, independently mutable selection-policy metadata.

    Returns:
        Policy separating native diagnostics from unchanged chemical checks.
    """
    return {
        "policy_version": "2026-11-24",
        "native_thr_role": NATIVE_THR_ROLE,
        "geometry_eligible_tiers": ["A", "B"],
        "native_thresholds_changed": False,
        "native_pass_and_joint_pass": "Historical diagnostic flags; no selection veto",
        "chemical_integrity_and_protonation_checks": "Unchanged",
    }
