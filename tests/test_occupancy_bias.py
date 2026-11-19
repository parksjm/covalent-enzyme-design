"""Synthetic checks for analytic bias and native evaluator unit contracts."""

from dataclasses import replace

import numpy as np
import pytest

from nyl12_covalent.occupancy_bias import (
    HARTREE_TO_KCAL_MOL,
    ActiveBiasCallback,
    OccupancyBiasConfig,
    evaluate_bias,
    evaluate_biased,
)


@pytest.fixture
def xyz() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.8, 0.2, 0.0],
            [1.4, 0.5, 0.1],
            [0.0, 2.0, 0.0],
            [0.0, 4.0, 0.2],
        ]
    )


@pytest.fixture
def config() -> OccupancyBiasConfig:
    return OccupancyBiasConfig(0, 1, (2, 3), (4, 5), 0.3, 1.7, 10.0, 4.0)


def test_known_energy(xyz: np.ndarray, config: OccupancyBiasConfig) -> None:
    initial = evaluate_bias(xyz, config)
    result = evaluate_bias(xyz, replace(config, q_target=initial.q - 0.2, k_r=0.0))
    assert result.energy == pytest.approx(0.2)
    assert result.q == pytest.approx(result.cn_acceptor - result.cn_donor)


def test_analytic_gradient(xyz: np.ndarray, config: OccupancyBiasConfig) -> None:
    result = evaluate_bias(xyz, config)
    numerical = np.empty_like(xyz)
    for index in np.ndindex(xyz.shape):
        plus, minus = xyz.copy(), xyz.copy()
        plus[index] += 1e-6
        minus[index] -= 1e-6
        numerical[index] = (
            evaluate_bias(plus, config).energy - evaluate_bias(minus, config).energy
        ) / 2e-6
    np.testing.assert_allclose(result.gradient, numerical, atol=2e-8)


def test_hydrogen_permutation(xyz: np.ndarray, config: OccupancyBiasConfig) -> None:
    first = evaluate_bias(xyz, config)
    permuted = xyz.copy()
    permuted[[2, 3]] = xyz[[3, 2]]
    second = evaluate_bias(permuted, config)
    assert second.energy == pytest.approx(first.energy)
    np.testing.assert_allclose(second.gradient[[3, 2]], first.gradient[[2, 3]])


def test_physical_units_and_active_callback(
    xyz: np.ndarray,
    config: OccupancyBiasConfig,
) -> None:
    seen = []

    def physical(positions: np.ndarray) -> tuple[float, np.ndarray]:
        seen.append(positions.copy())
        return float(np.sum(positions**2) / HARTREE_TO_KCAL_MOL), 2 * positions

    active = np.array([1, 2, 4])
    callback = ActiveBiasCallback(physical, xyz, active, config)
    positions = xyz[active].ravel() + 0.1
    energy, gradient = callback(positions)
    full = callback.coordinates(positions)
    assessment = evaluate_bias(full, config)
    assert energy == pytest.approx(np.sum(full**2) + assessment.energy)
    np.testing.assert_allclose(
        gradient, (2 * full + assessment.gradient)[active].ravel()
    )
    np.testing.assert_array_equal(seen[-1][[0, 3, 5]], xyz[[0, 3, 5]])
    assert callback.last_evaluation.physical_energy == pytest.approx(np.sum(full**2))
    numerical = np.empty_like(positions)
    for index in range(len(positions)):
        plus, minus = positions.copy(), positions.copy()
        plus[index] += 1e-6
        minus[index] -= 1e-6
        numerical[index] = (callback(plus)[0] - callback(minus)[0]) / 2e-6
    np.testing.assert_allclose(gradient, numerical, atol=2e-8)


@pytest.mark.parametrize(
    "changes",
    [
        {"k_q": -1},
        {"k_r": np.nan},
        {"r_target": 0},
        {"attack_atoms": (4, 4)},
        {"attack_atoms": (4, 8)},
        {"attack_atoms": (4.0, 5.0)},
        {"donor": 1},
        {"hydrogens": ()},
    ],
)
def test_invalid_config(
    xyz: np.ndarray, config: OccupancyBiasConfig, changes: dict
) -> None:
    with pytest.raises(ValueError):
        evaluate_bias(xyz, replace(config, **changes))


@pytest.mark.parametrize("active", [[], [1, 1], [-1], [7], [1.5], [[1]]])
def test_invalid_active(
    xyz: np.ndarray, config: OccupancyBiasConfig, active: list
) -> None:
    with pytest.raises(ValueError):
        ActiveBiasCallback(lambda x: (0.0, x), xyz, active, config)


def test_invalid_coordinates_and_physical(
    xyz: np.ndarray,
    config: OccupancyBiasConfig,
) -> None:
    zero_distance = xyz.copy()
    zero_distance[5] = zero_distance[4]
    with pytest.raises(ValueError):
        evaluate_bias(zero_distance, config)
    with pytest.raises(ValueError):
        evaluate_bias(np.full_like(xyz, np.nan), config)
    for energy, gradient in [(np.nan, xyz), (0.0, xyz[:2]), (0.0, xyz * np.nan)]:
        with pytest.raises(ValueError):
            evaluate_biased(xyz, lambda x: (energy, gradient), config)
    callback = ActiveBiasCallback(lambda x: (0.0, x), xyz, [1], config)
    for active in [np.zeros(4), np.full(3, np.nan), np.zeros((1, 3))]:
        with pytest.raises(ValueError):
            callback(active)
    assert callback.last_evaluation is None
