"""Regression checks for numerically stable relay path histograms."""

import numpy as np
import pytest

from nyl12_covalent.relay_umbrella_projection import path_histogram


def test_stable_accumulation_preserves_tiny_mass_and_sample_counts() -> None:
    result = path_histogram(
        np.array([0.1, 1.1, 4.0]),
        np.array([1.0, 1e-20, 2e-20]),
        np.arange(5.0),
        300.0,
        (0.4, 0.6),
        (2.4, 2.6),
    )

    assert result["probability_mass"] == [1.0, 1e-20, 0.0, 2e-20]
    assert result["raw_sample_count"] == [1, 1, 0, 1]


@pytest.mark.parametrize(
    "weights",
    [
        np.array([-0.1, 1.1]),
        np.array([np.nan, 1.0]),
        np.array([0.0, 0.0]),
        np.array([1.0]),
        np.ones((2, 1)),
    ],
)
def test_invalid_weights_fail_explicitly(weights: np.ndarray) -> None:
    with pytest.raises(ValueError, match="weights"):
        path_histogram(
            np.array([0.1, 0.9]),
            weights,
            np.linspace(0, 1, 6),
            300.0,
            (0.0, 0.2),
            (0.8, 1.0),
        )


def test_nonpositive_temperature_fails_explicitly() -> None:
    with pytest.raises(ValueError, match="Temperature"):
        path_histogram(
            np.array([0.1, 0.9]),
            np.array([0.5, 0.5]),
            np.linspace(0, 1, 6),
            0.0,
            (0.0, 0.2),
            (0.8, 1.0),
        )


def test_barrier_ignores_unsampled_outer_tails() -> None:
    centers = np.arange(0.15, 0.9, 0.1)
    weights = np.array([0.2, 0.4, 0.6, 0.3, 0.1, 0.02, 0.05, 0.01])
    result = path_histogram(
        centers,
        weights,
        np.arange(0.0, 1.1, 0.1),
        300.0,
        (0.2, 0.5),
        (0.6, 0.9),
    )

    pmf = np.array(
        [np.nan if value is None else value for value in result["pmf_kcal_mol"]]
    )
    assert result["status"] == "supported"
    assert np.isnan(pmf[[0, -1]]).all()
    assert np.isfinite(result["minimum_to_barrier_kcal_mol"])
    assert result["minimum_to_barrier_kcal_mol"] == pytest.approx(pmf[3:8].max())
