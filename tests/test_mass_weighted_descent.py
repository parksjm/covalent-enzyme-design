"""Tests for the pure mass-weighted descent pilot."""

import numpy as np
import pytest

from nyl12_covalent.mass_weighted_descent import mass_weighted_descent_step


def test_harmonic_1d_step_lowers_energy() -> None:
    def harmonic(x: np.ndarray) -> tuple[float, np.ndarray]:
        return 0.5 * float(np.sum(x**2)), x

    x0 = np.array([[1.0, 0.0, 0.0]])
    result = mass_weighted_descent_step(x0, np.array([1.0]), harmonic, 0.2)

    assert result.accepted
    assert result.energy < harmonic(x0)[0]
    assert result.coordinates[0, 0] == pytest.approx(0.8)


def test_harmonic_2d_path_follows_downhill_gradient() -> None:
    stiffness = np.array([1.0, 4.0, 0.0])

    def harmonic(x: np.ndarray) -> tuple[float, np.ndarray]:
        return 0.5 * float(np.sum(stiffness * x**2)), stiffness * x

    x0 = np.array([[1.0, 1.0, 0.0]])
    result = mass_weighted_descent_step(x0, np.ones(1), harmonic, 0.1)

    assert result.accepted
    assert result.coordinates[0, 0] < x0[0, 0]
    assert result.coordinates[0, 1] < x0[0, 1]
    assert result.energy < harmonic(x0)[0]


def test_mass_weighting_scales_cartesian_displacement_by_inverse_mass() -> None:
    gradient = np.array([[1.0, 1.0, 0.0]])

    def linear(x: np.ndarray) -> tuple[float, np.ndarray]:
        return float(np.sum(gradient * x)), gradient

    x0 = np.zeros((1, 3))
    result = mass_weighted_descent_step(x0, np.array([4.0]), linear, 0.2)
    expected = -0.2 / np.sqrt(2.0) / 2.0

    assert result.accepted
    assert result.coordinates[0, 0] == pytest.approx(expected)
    assert result.coordinates[0, 1] == pytest.approx(expected)


def test_energy_increase_rejects_and_halves_step() -> None:
    def harmonic(x: np.ndarray) -> tuple[float, np.ndarray]:
        return 0.5 * float(np.sum(x**2)), x

    x0 = np.array([[1.0, 0.0, 0.0]])
    result = mass_weighted_descent_step(
        x0,
        np.ones(1),
        harmonic,
        0.2,
        initial_direction_q=np.array([[1.0, 0.0, 0.0]]),
    )

    assert not result.accepted
    np.testing.assert_array_equal(result.coordinates, x0)
    assert result.energy == harmonic(x0)[0]
    assert result.step_size == pytest.approx(0.1)


def test_fixed_atoms_remain_unchanged() -> None:
    def harmonic(x: np.ndarray) -> tuple[float, np.ndarray]:
        return 0.5 * float(np.sum(x**2)), x

    x0 = np.array([[1.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
    result = mass_weighted_descent_step(
        x0, np.ones(2), harmonic, 0.1, fixed_atoms=np.array([False, True])
    )

    assert result.accepted
    np.testing.assert_array_equal(result.coordinates[1], x0[1])


@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_both_branches_can_be_selected_at_exact_saddle(sign: float) -> None:
    def saddle(x: np.ndarray) -> tuple[float, np.ndarray]:
        return -0.5 * float(x[0, 0] ** 2), np.array([[-x[0, 0], 0.0, 0.0]])

    x0 = np.zeros((1, 3))
    result = mass_weighted_descent_step(
        x0,
        np.ones(1),
        saddle,
        0.1,
        initial_direction_q=np.array([[sign, 0.0, 0.0]]),
    )

    assert result.accepted
    assert np.sign(result.coordinates[0, 0]) == sign


def test_nonfinite_gradient_fails_explicitly() -> None:
    def invalid(x: np.ndarray) -> tuple[float, np.ndarray]:
        return 0.0, np.full_like(x, np.nan)

    with pytest.raises(ValueError, match="finite"):
        mass_weighted_descent_step(np.zeros((1, 3)), np.ones(1), invalid, 0.1)
