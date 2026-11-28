"""Native trajectory time-grid validation tests."""

import numpy as np
import pytest
from pathlib import Path

from nyl12_covalent.relay_umbrella_diagnostics import (
    _read_dumpave,
    _step_origin_from_mdin,
    _validate_trajectory_times,
)


def test_identical_unsampled_startup_evaluations_are_coalesced(tmp_path: Path) -> None:
    """Fresh Amber stages can print the same initial restraint evaluation twice."""
    path = tmp_path / "umbrella_coordinate.dat"
    path.write_text("0 1.802 1.755\n0 1.802 1.755\n10 1.820 1.773\n")
    assert _read_dumpave(path) == {0: (1.802, 1.755), 10: (1.820, 1.773)}


@pytest.mark.parametrize(
    "records",
    ["0 1.802 1.755\n0 1.803 1.755\n", "10 1.802 1.755\n10 1.802 1.755\n"],
)
def test_conflicting_or_sampled_duplicate_cv_records_fail(
    tmp_path: Path, records: str
) -> None:
    """Only identical unsampled startup records are accepted as duplicates."""
    path = tmp_path / "umbrella_coordinate.dat"
    path.write_text(records)
    with pytest.raises(ValueError, match="Duplicate/nonfinite"):
        _read_dumpave(path)


def native_times() -> np.ndarray:
    """Return float32 Amber-like absolute times at 100-step intervals."""
    steps = np.arange(100, 10_001, 100, dtype=np.float64)
    return np.asarray(6.0 + steps * 0.0005, dtype=np.float32)


def test_float32_serialization_error_is_bounded_and_steps_are_exact() -> None:
    times = native_times()
    assert times.dtype == np.float32
    assert float(times[0]) == pytest.approx(6.05)

    result = _validate_trajectory_times(times, 12_000, 10_000, 100)

    assert result == list(range(100, 10_001, 100))


def test_float64_times_use_only_small_arithmetic_allowance() -> None:
    times = np.asarray(6.0 + np.arange(1, 101) * 0.05, dtype=np.float64)

    assert _validate_trajectory_times(times, 12_000, 10_000, 100) == list(
        range(100, 10_001, 100)
    )


def test_fresh_restart_uses_zero_clock_even_when_seed_is_old() -> None:
    mdin = "&cntrl irest=0, ntx=1, nstlim=2000, dt=0.0005, ntwx=100, /"
    times = np.asarray(np.arange(1, 21, dtype=np.float64) * 0.05, dtype=np.float32)

    origin = _step_origin_from_mdin(mdin, 36.0)

    assert origin == 0
    assert _validate_trajectory_times(times, origin, 2_000, 100) == list(
        range(100, 2_001, 100)
    )


def test_continuing_restart_uses_restart_clock() -> None:
    mdin = "&cntrl irest=1, ntx=5, nstlim=2000, dt=0.0005, ntwx=100, /"
    times = np.asarray(1.0 + np.arange(1, 21) * 0.05, dtype=np.float32)

    origin = _step_origin_from_mdin(mdin, 1.0)

    assert origin == 2_000
    assert _validate_trajectory_times(times, origin, 2_000, 100) == list(
        range(100, 2_001, 100)
    )


@pytest.mark.parametrize(
    ("mdin", "restart_time_ps", "message"),
    [
        ("&cntrl ntx=1, /", 0.0, "declare irest"),
        ("&cntrl irest=2, /", 0.0, "Unsupported Amber irest"),
        ("&cntrl irest=1.0, /", 0.0, "must be integer"),
        ("&cntrl irest=0, irest=1, /", 0.0, "declare irest exactly once"),
        ("&cntrl irest=1, /", float("nan"), "Restart time"),
    ],
)
def test_clock_origin_rejects_invalid_or_ambiguous_input(
    mdin: str, restart_time_ps: float, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        _step_origin_from_mdin(mdin, restart_time_ps)


def test_time_outside_native_ulp_bound_is_rejected() -> None:
    times = native_times()
    times[0] = np.float32(times[0] + 0.00001)

    with pytest.raises(ValueError, match="off the integration grid"):
        _validate_trajectory_times(times, 12_000, 10_000, 100)


def test_wrong_sampled_cadence_is_rejected() -> None:
    times = native_times()
    times[1] = np.float32(times[1] + 0.0005)

    with pytest.raises(ValueError, match="off the integration grid"):
        _validate_trajectory_times(times, 12_000, 10_000, 100)


def test_nonmonotonic_times_are_rejected() -> None:
    times = native_times()
    times[1] = times[0]

    with pytest.raises(ValueError, match="strictly increasing"):
        _validate_trajectory_times(times, 12_000, 10_000, 100)


def test_frame_count_and_stride_must_match_completed_steps() -> None:
    with pytest.raises(ValueError, match="expected 100"):
        _validate_trajectory_times(native_times()[:-1], 12_000, 10_000, 100)
    with pytest.raises(ValueError, match="divisible"):
        _validate_trajectory_times(native_times(), 12_000, 10_001, 100)
    with pytest.raises(ValueError, match="divisible"):
        _validate_trajectory_times(native_times(), 12_000, 0, 100)
    with pytest.raises(ValueError, match="divisible"):
        _validate_trajectory_times(native_times(), 12_000, 10_000, 0)


def test_times_must_be_finite_one_dimensional_floats() -> None:
    with pytest.raises(ValueError, match="one-dimensional float"):
        _validate_trajectory_times(np.array([1, 2]), 0, 2, 1)
    with pytest.raises(ValueError, match="one-dimensional float"):
        _validate_trajectory_times(np.array([[0.05, 0.1]]), 0, 2, 1)
    times = native_times()
    times[4] = np.nan
    with pytest.raises(ValueError, match="finite"):
        _validate_trajectory_times(times, 12_000, 10_000, 100)
