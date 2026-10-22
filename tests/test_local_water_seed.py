"""Tests for deterministic local water placement."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.local_water_seed import seed_water


def test_preserve_original_oxygen_allows_finite_close_contacts() -> None:
    oxygen = np.array([0.0, 0.0, 0.0])
    nitrogen, og1 = np.array([-2.8, 0.0, 0.0]), np.array([0.0, 2.8, 0.0])
    environment = np.array([nitrogen, og1, [1.67, 0.0, 0.0]])
    water, report = seed_water(
        oxygen, nitrogen, og1, environment, np.array([7, 8, 6]), preserve_oxygen=True
    )
    assert np.array_equal(water[0], oxygen)
    assert report["oxygen_displacement_A"] == 0
    assert report["minimum_heavy_distance_A"] == pytest.approx(1.67)


def test_preserve_oxygen_rejects_exact_atom_coincidence() -> None:
    oxygen = np.zeros(3)
    nitrogen, og1 = np.array([-2.8, 0.0, 0.0]), np.array([0.0, 2.8, 0.0])
    with pytest.raises(ValueError, match="overlap"):
        seed_water(
            oxygen,
            nitrogen,
            og1,
            np.array([nitrogen, og1, oxygen]),
            np.array([7, 8, 6]),
            preserve_oxygen=True,
        )


def _bridge_geometry() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a symmetric pair of anchors and their ideal bridge oxygen."""
    nitrogen = np.array([0.0, 0.0, 0.0])
    og1 = np.array([5.6, 0.0, 0.0])
    oxygen = np.array([2.8, 0.0, 0.0])
    return oxygen, nitrogen, og1


def test_seed_water_builds_correct_bridge_and_water_geometry() -> None:
    oxygen, nitrogen, og1 = _bridge_geometry()
    environment = np.array([[20.0, 20.0, 20.0]])

    water, _ = seed_water(oxygen, nitrogen, og1, environment, np.array([6]))

    assert water.shape == (3, 3)
    assert np.all(np.isfinite(water))
    for anchor in (nitrogen, og1):
        assert 2.6 <= np.linalg.norm(water[0] - anchor) <= 3.2
    assert np.linalg.norm(water[1] - water[0]) == pytest.approx(0.9572, abs=1e-6)
    assert np.linalg.norm(water[2] - water[0]) == pytest.approx(0.9572, abs=1e-6)
    first = water[1] - water[0]
    second = water[2] - water[0]
    cosine = np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
    angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    assert angle == pytest.approx(104.52, abs=1e-5)
    assert np.linalg.norm(water[0] - oxygen) <= 1.5


def test_seed_water_orientation_is_deterministic() -> None:
    oxygen, nitrogen, og1 = _bridge_geometry()
    environment = np.array([[20.0, 20.0, 20.0]])

    elements = np.array([6])
    first, _ = seed_water(oxygen, nitrogen, og1, environment, elements)
    second, _ = seed_water(oxygen, nitrogen, og1, environment, elements)

    np.testing.assert_array_equal(first, second)


def test_seed_water_rejects_missing_bridge_candidate() -> None:
    oxygen, nitrogen, _ = _bridge_geometry()
    distant_og1 = np.array([20.0, 0.0, 0.0])

    with pytest.raises(ValueError, match="bridge ranges"):
        seed_water(
            oxygen, nitrogen, distant_og1, np.empty((0, 3)), np.array([], dtype=int)
        )


def test_seed_water_rejects_catastrophic_heavy_atom_overlap() -> None:
    oxygen, nitrogen, og1 = _bridge_geometry()
    environment = np.array([oxygen.copy()])

    with pytest.raises(ValueError, match="overlap|distance|candidate"):
        seed_water(oxygen, nitrogen, og1, environment, np.array([6]))


@pytest.mark.parametrize(
    ("environment", "elements"),
    [
        (np.zeros((2, 2)), np.array([6, 6])),
        (np.zeros((2, 3)), np.array([6])),
        (np.array([[np.nan, 0.0, 0.0]]), np.array([6])),
    ],
)
def test_seed_water_rejects_invalid_environment(
    environment: np.ndarray, elements: tuple[str, ...]
) -> None:
    oxygen, nitrogen, og1 = _bridge_geometry()

    with pytest.raises(ValueError):
        seed_water(oxygen, nitrogen, og1, environment, elements)


def test_seed_water_rejects_unsupported_environment_element() -> None:
    oxygen, nitrogen, og1 = _bridge_geometry()

    with pytest.raises(KeyError):
        seed_water(
            oxygen,
            nitrogen,
            og1,
            np.array([[20.0, 20.0, 20.0]]),
            np.array([999]),
        )


def test_seed_water_does_not_modify_inputs() -> None:
    oxygen, nitrogen, og1 = _bridge_geometry()
    environment = np.array([[20.0, 20.0, 20.0]])
    original = tuple(value.copy() for value in (oxygen, nitrogen, og1, environment))
    elements = np.array([6])

    seed_water(oxygen, nitrogen, og1, environment, elements)

    for value, expected in zip((oxygen, nitrogen, og1, environment), original):
        np.testing.assert_array_equal(value, expected)
    np.testing.assert_array_equal(elements, [6])
