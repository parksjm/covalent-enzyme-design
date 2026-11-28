"""Recoverable SCC warnings remain distinct from invalid returned gradients."""

import pytest

from nyl12_covalent.relay_scc_guard import scc_guard_reason
from nyl12_covalent.relay_ts_guarded import invalid_trial

WARNING = "QMMM SCC-DFTB: SCC-DFTB FOR STEP 1 DID NOT CONVERGE AFTER 70 cycles.\n"
RECOVERED = (
    WARNING
    + "QMMM SCC-DFTB: 96 -123.0 0.0 0.0\nQMMM SCC-DFTB: SCC-DFTB for step 1 converged in 96 cycles.\n"
)


def test_recovered_warning_does_not_reject_valid_evaluation() -> None:
    """Native70-to96cycle recovery is accepted with identical physical tolerances."""
    assert scc_guard_reason(RECOVERED, final=True) is None
    assert invalid_trial(RECOVERED, -718.0, 25.0, final=True) is None
    assert scc_guard_reason(WARNING + "QMMM SCC-DFTB: 82 -123.0 0.0 0.0\n") is None
    assert (
        scc_guard_reason(RECOVERED + "Energy calculation finished, energy: -718.0")
        is None
    )


@pytest.mark.parametrize(
    "tail",
    [
        "QMMM SCC-DFTB: 501 -123.0 1.0 1.0\n",
        "QMMM SCC-DFTB: SCC-DFTB FOR STEP 1 DID NOT CONVERGE AFTER 500 cycles.\n",
        "QMMM SCC-DFTB: SCC-DFTB for step 1 converged in 501 cycles.\n",
    ],
)
def test_long_oscillation_exhausts_budget(tail: str) -> None:
    """A500iteration budget caps active oscillation and over-budget recovery."""
    assert "500-iteration" in scc_guard_reason(WARNING + tail)


def test_fatal_and_final_unresolved_fail_immediately() -> None:
    """Pending recovery is allowed only while the electronic trial remains active."""
    assert scc_guard_reason(WARNING, final=True) is not None
    assert scc_guard_reason("QM calculation failed") is not None
    assert scc_guard_reason("increasing the electronic temperature") is not None
    assert (
        scc_guard_reason(WARNING + "Energy calculation finished, energy: -718.0")
        is not None
    )


def test_later_same_step_label_cannot_mask_failed_returned_gradient() -> None:
    """Repeated native step1labels cannot heal an earlier invalid force evaluation."""
    text = WARNING + "Energy calculation finished, energy: -718.0\n" + RECOVERED
    assert "unresolved SCC electronic gradient" in scc_guard_reason(text)
    with pytest.raises(ValueError):
        scc_guard_reason("", maximum_iterations=0)
