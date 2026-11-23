"""Tests for native DFTB SCC history classification."""

import pytest

from nyl12_covalent.dftb_scc_validation import audit_dftb_scc


def trial_failed(cycles: int, step: int = 1) -> str:
    """Return one native recoverable trial warning."""
    return (
        " QMMM SCC-DFTB: SCC-DFTB FOR STEP "
        f"{step} DID NOT CONVERGE AFTER {cycles} cycles."
    )


def trial_converged(cycles: int, step: int = 1) -> str:
    """Return one native SCC convergence summary."""
    return f" QMMM SCC-DFTB: SCC-DFTB for step {step} converged in {cycles} cycles."


def test_exact_same_step_recoveries_record_each_trial() -> None:
    """The exact retry warning is recoverable only after a same-step success."""
    result = audit_dftb_scc(
        "\n".join(
            [
                "QMMM SCC-DFTB:      Telec :     0.000K",
                trial_failed(70),
                trial_converged(14),
                trial_failed(140),
                trial_converged(22),
                trial_failed(210),
                trial_converged(31),
            ]
        )
    )

    assert result["passed"]
    assert result["history_observed"]
    assert result["verified_scc_history"]
    assert result["recovered_trial_count"] == 3
    assert result["recovered_trials"] == [
        {"step": 1, "failed_cycles": 70, "recovered_cycles": 14},
        {"step": 1, "failed_cycles": 140, "recovered_cycles": 22},
        {"step": 1, "failed_cycles": 210, "recovered_cycles": 31},
    ]


@pytest.mark.parametrize(
    "text",
    [
        "QMMM SCC-DFTB: Convergence could not be achieved in this step.",
        "QMMM SCC-DFTB: SCC Convergence failure - ITRMAX exceeded.",
        "QMMM SCC-DFTB: **WARNING** SCC-DFTB FOR STEP 1 DID NOT CONVERGE AFTER 70 cycles.",
        "QMMM SCF failed to converge.",
        "QMMM SCC-DFTB: forces will not be accurate.",
        "QM calculation failed.",
        "Final energy = NaN",
    ],
)
def test_fatal_or_nonwhitelisted_failure_is_rejected(text: str) -> None:
    """Only the exact native inner-trial warning can enter retry recovery."""
    result = audit_dftb_scc(text)

    assert not result["passed"]
    assert result["fatal_markers"]


def test_failure_after_recovery_remains_fatal() -> None:
    """A later genuine convergence failure is not erased by earlier recovery."""
    result = audit_dftb_scc(
        "\n".join(
            [
                trial_failed(70),
                trial_converged(12),
                "QMMM SCC-DFTB: Convergence could not be achieved in this step.",
            ]
        )
    )

    assert not result["passed"]
    assert result["recovered_trial_count"] == 1
    assert result["fatal_markers"]


def test_unrecovered_and_different_step_warnings_are_rejected() -> None:
    """Warnings remain pending until their own step has a later success."""
    dangling = audit_dftb_scc(trial_failed(70))
    other_step = audit_dftb_scc(
        "\n".join([trial_failed(70, step=2), trial_converged(9)])
    )

    assert not dangling["passed"]
    assert dangling["unresolved_trials"] == [{"step": 1, "failed_cycles": 70}]
    assert not other_step["passed"]
    assert other_step["unresolved_trials"] == [{"step": 2, "failed_cycles": 70}]


@pytest.mark.parametrize(
    "text",
    [
        "QMMM SCC-DFTB:      Telec :   300.000K",
        "QMMM SCC-DFTB: Increasing the electronic temperature of step 1 to 300 K.",
        "QMMM SCC-DFTB: The energy for this step was calculated with a higher temperature.",
        "QMMM SCC-DFTB: Electronic temperature restored to 0 K.",
    ],
)
def test_nonzero_or_changed_electronic_temperature_is_rejected(text: str) -> None:
    """Temperature changes remain outside the accepted zero-K history."""
    result = audit_dftb_scc(text)

    assert not result["passed"]
    assert result["temperature_changes"]


def test_zero_temperature_and_marker_free_history_are_distinguished() -> None:
    """Zero-K evidence passes, while absent history is never called verified."""
    zero = audit_dftb_scc("QMMM SCC-DFTB:      Telec :     0.000K")
    marker_free = audit_dftb_scc(
        "Optimisation converged\nworking directory: /tmp/banana"
    )

    assert zero["passed"]
    assert not zero["verified_scc_history"]
    assert marker_free["passed"]
    assert not marker_free["history_observed"]
    assert not marker_free["verified_scc_history"]
    assert marker_free["fatal_markers"] == []


def test_nan_word_boundary_rejects_nan_but_not_directory_substring() -> None:
    """NaN is matched as a token rather than text inside an ordinary path."""
    assert not audit_dftb_scc("energy=NaN")["passed"]
    assert audit_dftb_scc("working directory: /tmp/banana")["passed"]
