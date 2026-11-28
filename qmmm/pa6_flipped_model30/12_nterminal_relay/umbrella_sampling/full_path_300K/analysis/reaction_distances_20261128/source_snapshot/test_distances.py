"""Regression for weighted empirical distributions with one dominant sample."""

import numpy as np
from run_distances import summarize


def test_dominant_sample_preserves_empirical_quantiles() -> None:
    """A dominant observation controls quantiles without interpolated geometry."""
    result = summarize(np.array([5.0, 1.0, 3.0]), np.array([0.97, 0.02, 0.01]))
    assert result["p10_A"] == result["median_A"] == result["p90_A"] == 5.0
    assert np.isclose(result["mean_A"], 4.9)
    uniform = summarize(np.array([3.0, 1.0, 2.0]), np.ones(3))
    assert (uniform["p10_A"], uniform["median_A"], uniform["p90_A"]) == (1.0, 2.0, 3.0)
