"""Validate label invariance and the custom occupancy derivative."""

import numpy as np
import pytest

from nyl12_covalent.proton_coordination import (
    coordination_transfer,
    smooth_coordination,
)


def test_invariant_under_hydrogen_exchange() -> None:
    xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.7, 0.0, 0.0],
            [0.98, 0.0, 0.0],
            [3.2, 0.8, 0.0],
            [3.2, -0.8, 0.0],
            [0.2, 2.0, 0.0],
        ]
    )
    first = coordination_transfer(xyz, 0, 1, [2, 3, 4, 5])
    exchanged = xyz.copy()
    exchanged[[2, 5]] = exchanged[[5, 2]]
    second = coordination_transfer(exchanged, 0, 1, [2, 3, 4, 5])
    assert first.value == pytest.approx(second.value)
    assert first.donor == pytest.approx(second.donor)
    assert first.acceptor == pytest.approx(second.acceptor)
    expected = first.gradient.copy()
    expected[[2, 5]] = expected[[5, 2]]
    np.testing.assert_allclose(second.gradient, expected, atol=1e-14)


def test_custom_analytic_gradient() -> None:
    xyz = np.random.default_rng(42).normal(size=(7, 3))
    result = coordination_transfer(xyz, 0, 1, [2, 3, 4, 5, 6])
    numerical = np.zeros_like(xyz)
    step = 1e-6
    for atom, axis in np.ndindex(xyz.shape):
        plus, minus = xyz.copy(), xyz.copy()
        plus[atom, axis] += step
        minus[atom, axis] -= step
        numerical[atom, axis] = (
            coordination_transfer(plus, 0, 1, [2, 3, 4, 5, 6]).value
            - coordination_transfer(minus, 0, 1, [2, 3, 4, 5, 6]).value
        ) / (2 * step)
    np.testing.assert_allclose(result.gradient, numerical, atol=2e-9)
    np.testing.assert_allclose(result.gradient.sum(axis=0), 0, atol=1e-14)


@pytest.mark.parametrize(
    "center, hydrogens, radius",
    [
        (0, [], 1.25),
        (0, [1, 1], 1.25),
        (0, [0], 1.25),
        (3, [1], 1.25),
        (0, [3], 1.25),
        (0, [-1], 1.25),
        (0, [1.2], 1.25),
        (0, [1], 0),
        (0, [1], np.nan),
    ],
)
def test_invalid_inputs(center: int, hydrogens: list, radius: float) -> None:
    with pytest.raises(ValueError):
        smooth_coordination(np.zeros((3, 3)), center, hydrogens, radius)


def test_invalid_geometry_and_identical_sites() -> None:
    with pytest.raises(ValueError):
        smooth_coordination(np.zeros(3), 0, [1])
    with pytest.raises(ValueError):
        smooth_coordination(np.full((3, 3), np.nan), 0, [1])
    with pytest.raises(ValueError):
        coordination_transfer(np.zeros((3, 3)), 0, 0, [1])
