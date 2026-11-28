"""Focused tests for modular E-S QM/MM endpoint measurements."""

from typing import Any

import pytest

from scripts.queue_modular_es_qmmm_after_mm import (
    mm_acceptance_decision,
    native_gate,
    validate_override_options,
)


def test_native_gate_requires_every_predeclared_limit() -> None:
    thr = {
        "heavy_atom_displacements_angstrom": {"N": 1.0, "OG1": 0.7},
        "chi1_delta_es_minus_reference_degrees": 19.0,
        "psi_delta_es_minus_reference_degrees": -34.0,
    }
    report = {"states": {"continuation_endpoint": {"thr233_held_out": thr}}}
    assert native_gate(report)["accepted"] is True
    thr["heavy_atom_displacements_angstrom"]["OG1"] = 0.8
    assert native_gate(report)["accepted"] is False


def _acceptance_inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    return (
        {"accepted": True},
        {"tier": "B"},
        {
            "accepted": False,
            "checks": {
                "N_displacement_le_1.25_A": True,
                "OG1_displacement_le_0.75_A": False,
                "abs_chi1_delta_le_20_deg": True,
                "abs_psi_delta_le_35_deg": True,
            },
            "metrics": {"heavy_atom_displacements_angstrom": {"OG1": 0.782}},
        },
    )


def test_native_rejection_is_diagnostic_without_override() -> None:
    decision = mm_acceptance_decision(*_acceptance_inputs(), None, None)

    assert decision["standard_accepted_for_qmmm"] is True
    assert decision["legacy_native_filtered_accepted"] is False
    assert decision["selection_policy"]["native_thr_role"] == "diagnostic_only"
    assert decision["override"] == {
        "requested": False,
        "accepted": False,
        "reason": None,
        "og1_displacement_max_angstrom": None,
        "failed_conditions": [],
    }
    assert decision["accepted_for_qmmm"] is True


def test_authorized_og1_near_pass_accepts_and_preserves_standard_rejection() -> None:
    decision = mm_acceptance_decision(
        *_acceptance_inputs(), 0.80, "Reviewed converged Tier-B near pass"
    )

    assert decision["standard_accepted_for_qmmm"] is True
    assert decision["legacy_native_filtered_accepted"] is False
    assert decision["override"]["accepted"] is True
    assert decision["override"]["reason"] == "Reviewed converged Tier-B near pass"
    assert decision["override"]["og1_displacement_max_angstrom"] == 0.80
    assert decision["accepted_for_qmmm"] is True


@pytest.mark.parametrize(
    ("maximum", "reason"), [(0.80, None), (None, "explicit reason")]
)
def test_override_options_must_be_supplied_as_a_pair(
    maximum: float | None, reason: str | None
) -> None:
    with pytest.raises(ValueError, match="must be supplied together"):
        validate_override_options(maximum, reason)


def test_other_native_failures_also_remain_diagnostic() -> None:
    convergence, geometry, native = _acceptance_inputs()
    native["checks"]["abs_psi_delta_le_35_deg"] = False

    decision = mm_acceptance_decision(
        convergence, geometry, native, 0.80, "Reviewed near pass"
    )

    assert decision["override"]["accepted"] is False
    assert "abs_psi_delta_le_35_deg" in decision["override"]["failed_conditions"]
    assert decision["accepted_for_qmmm"] is True


@pytest.mark.parametrize("failure", ["convergence", "geometry"])
def test_diagnostic_only_policy_keeps_numerical_and_chemical_gates(
    failure: str,
) -> None:
    """Removing a native veto cannot admit a failed optimization or bad pose."""
    convergence, geometry, native = _acceptance_inputs()
    if failure == "convergence":
        convergence["accepted"] = False
    else:
        geometry["tier"] = "C"
    decision = mm_acceptance_decision(convergence, geometry, native, None, None)
    assert decision["accepted_for_qmmm"] is False
