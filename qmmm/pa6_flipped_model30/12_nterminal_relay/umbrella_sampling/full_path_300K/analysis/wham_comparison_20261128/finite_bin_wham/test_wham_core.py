"""Analytic qualification of the existing-solver finite-bin WHAM adapter."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import logsumexp

from wham_core import projected_pmf, quantize, solve_wham


def test_single_biased_state() -> None:
    """One-state WHAM must exactly remove the known cell bias."""
    u = np.array([[0.0, 0.7, 1.2]])
    inverse = np.repeat(np.arange(3), [20, 30, 40])
    got = solve_wham(u, inverse, np.array([90]))
    expected = np.array([20, 30, 40]) * np.exp(u[0])
    expected /= expected.sum()
    np.testing.assert_allclose(np.exp(got["cell_logmass"]), expected, atol=1e-12)
    assert got["fixed_point_residual_kBT"] < 1e-10


def fixture() -> tuple:
    """Construct exact integer histograms for known discrete unbiased masses."""
    p = np.array([0.2, 0.3, 0.5])
    biased = np.array([[0.5, 0.3, 0.2], [0.1, 0.3, 0.6]])
    u = -np.log(biased / p)
    inverse = np.r_[
        np.repeat(np.arange(3), [50, 30, 20]),
        np.repeat(np.arange(3), [10, 30, 60]),
    ]
    return p, u, inverse, np.array([100, 100])


def test_multistate_and_gauge() -> None:
    """Known masses and state-specific energy-gauge invariance are recovered."""
    p, u, inverse, counts = fixture()
    got = solve_wham(u, inverse, counts)
    np.testing.assert_allclose(np.exp(got["cell_logmass"]), p, atol=1e-10)
    shifted = solve_wham(u + np.array([[5.0], [-3.0]]), inverse, counts)
    np.testing.assert_allclose(got["logweights"], shifted["logweights"], atol=1e-9)
    assert abs(logsumexp(got["logweights"])) < 1e-12


def test_identical_bias_split() -> None:
    """Splitting a state's sample counts does not change the point estimate."""
    _, u, inverse, counts = fixture()
    original = solve_wham(u, inverse, counts)
    split = solve_wham(u[[0, 0, 1]], inverse, np.array([40, 60, 100]))
    np.testing.assert_allclose(original["logweights"], split["logweights"], atol=1e-9)


def test_equal_cell_weights_and_projection() -> None:
    """Original-frame marginal equals cell mass times empirical conditional bins."""
    _, u, inverse, counts = fixture()
    got = solve_wham(u, inverse, counts)
    frame_w = np.exp(got["logweights"])
    for cell in range(3):
        assert np.ptp(frame_w[inverse == cell]) == 0
    values = np.arange(len(inverse)) % 2 + 0.5
    edges = np.array([0, 1, 2])
    _, mass = projected_pmf(values, got["logweights"], edges, 0, 0.6)
    expected = np.zeros(2)
    for cell, p in enumerate(np.exp(got["cell_logmass"])):
        expected += (
            p
            * np.bincount(values[inverse == cell].astype(int), minlength=2)
            / (inverse == cell).sum()
        )
    np.testing.assert_allclose(mass, expected, atol=1e-12)
    assert abs(mass.sum() - 1) < 1e-12


def test_quantize_and_fixed_reference() -> None:
    """Cell centers respect negative coordinates and explicit half-cell shifts."""
    q = np.array([[-0.01, 0.01], [0.02, 0.03], [0.09, 0.08]])
    centers, inverse, occupancy = quantize(q, 0.1)
    np.testing.assert_allclose(centers, [[-0.05, 0.05], [0.05, 0.05]])
    np.testing.assert_array_equal(inverse, [0, 1, 1])
    np.testing.assert_array_equal(occupancy, [1, 2])
    shifted, _, _ = quantize(q, 0.1, 0.5)
    np.testing.assert_allclose(shifted, [[0, 0], [0.1, 0.1]])
    pmf, mass = projected_pmf(
        np.array([0.1, 0.2, 1.5]), np.log([0.25, 0.25, 0.5]), np.arange(4), 0, 1
    )
    np.testing.assert_allclose(pmf[:2], 0)
    assert np.isnan(pmf[2])
    np.testing.assert_allclose(mass, [0.5, 0.5, 0])


@pytest.mark.parametrize("width,shift", [(0, 0), (-1, 0), (np.nan, 0), (1, np.inf)])
def test_invalid_grid(width: float, shift: float) -> None:
    """Invalid scientific discretizations fail explicitly."""
    with pytest.raises(ValueError):
        quantize(np.ones((2, 3)), width, shift)


def test_invalid_inputs() -> None:
    """Unsupported states, empty cells and missing reference support fail."""
    with pytest.raises(ValueError):
        quantize(np.array([[np.nan]]), 1)
    with pytest.raises(ValueError):
        solve_wham(np.zeros((1, 1)), np.array([0]), np.array([0]))
    with pytest.raises(ValueError):
        solve_wham(np.zeros((1, 2)), np.array([0]), np.array([1]))
    with pytest.raises(ValueError):
        projected_pmf(np.array([1.5]), np.zeros(1), np.arange(3), 0, 1)
    with pytest.raises(ValueError):
        projected_pmf(np.array([3.5]), np.zeros(1), np.arange(3), 0, 1)
