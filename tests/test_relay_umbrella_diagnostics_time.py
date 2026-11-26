"""Native trajectory time-grid validation tests."""

import numpy as np
import pytest

from nyl12_covalent.relay_umbrella_diagnostics import (
    _validate_trajectory_times,
)


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
