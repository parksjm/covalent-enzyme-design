"""Strict electronic history and total-energy divergence guard checks."""

import pytest

from nyl12_covalent.relay_ts_guarded import invalid_trial


@pytest.mark.parametrize(
    "text",
    [
        "QM calculation failed",
        "increasing the electronic temperature",
        "Etot = NaN\nQM calculation failed",
    ],
)
def test_invalid_electronic_trials_abort(text: str) -> None:
    """Even recovered SCC retries cannot continue a strict repaired search."""
    assert invalid_trial(text, -718.0, 25.0) is not None


def test_energy_guard_uses_absolute_hartree_to_native_conversion() -> None:
    """The fixed declared rise guard detects runaway without geometric walls."""
    assert invalid_trial("Etot = -450000.0", -718.0, 25.0) is not None
    assert invalid_trial("Etot = -450560.0", -718.0, 25.0) is None
    assert (
        invalid_trial(
            "QMMM SCC-DFTB: SCC-DFTB for step 1 converged in 10 cycles.", -718.0, 25.0
        )
        is None
    )
