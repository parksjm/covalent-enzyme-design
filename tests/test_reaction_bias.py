"""Analytic-gradient, unit, permutation and input-contract checks."""

import json
from dataclasses import replace

import numpy as np
import pytest

from nyl12_covalent.reaction_bias import (
    HARTREE_TO_KCAL_MOL,
    CoordinationCV,
    CoordinationDifferenceCV,
    DistanceCV,
    HarmonicTarget,
    ReactionBiasCallback,
    evaluate_bias,
    evaluate_biased,
)


@pytest.fixture
def xyz() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.1, 0.0],
            [3.0, 1.0, 0.2],
            [0.8, 0.2, 0.0],
            [1.4, 0.5, 0.1],
            [2.8, 0.8, 0.3],
        ]
    )


@pytest.fixture
def terms() -> tuple[HarmonicTarget, ...]:
    return (
        HarmonicTarget(DistanceCV(0, 1), 1.5, 8.0),
        HarmonicTarget(DistanceCV(1, 2), 2.0, 9.0),
        HarmonicTarget(CoordinationDifferenceCV(0, 1, (3, 4, 5)), 0.5, 10.0),
        HarmonicTarget(CoordinationDifferenceCV(1, 2, (3, 4, 5)), -0.2, 7.0),
        HarmonicTarget(CoordinationCV(2, (3, 4, 5)), 1.0, 3.0),
    )


def test_gradient_energy_force_and_metadata(xyz: np.ndarray, terms: tuple) -> None:
    result = evaluate_bias(xyz, terms)
    assert len(result.values) == len(result.targets) == len(result.energies) == 5
    assert result.energy == pytest.approx(sum(result.energies))
    assert result.targets == tuple(term.target for term in terms)
    np.testing.assert_allclose(result.gradient.sum(axis=0), 0.0, atol=1e-14)
    numerical = np.empty_like(xyz)
    for index in np.ndindex(xyz.shape):
        plus, minus = xyz.copy(), xyz.copy()
        plus[index] += 1e-6
        minus[index] -= 1e-6
        numerical[index] = (
            evaluate_bias(plus, terms).energy - evaluate_bias(minus, terms).energy
        ) / 2e-6
    np.testing.assert_allclose(result.gradient, numerical, atol=3e-8)
    metadata = json.loads(json.dumps([term.to_metadata() for term in terms]))
    assert metadata[2]["cv"]["type"] == "CoordinationDifferenceCV"
    assert metadata[2]["cv"]["hydrogens"] == [3, 4, 5]
    known = HarmonicTarget(DistanceCV(0, 1), result.values[0] - 0.2, 10.0)
    assert evaluate_bias(xyz, [known]).energy == pytest.approx(0.2)


def test_hydrogen_permutation(xyz: np.ndarray, terms: tuple) -> None:
    first = evaluate_bias(xyz, terms)
    order = [0, 1, 2, 5, 3, 4]
    second = evaluate_bias(xyz[order], terms)
    assert second.energy == pytest.approx(first.energy)
    np.testing.assert_allclose(second.gradient, first.gradient[order], atol=1e-14)
    # All 54 real QM hydrogens are accepted without an identity restriction.
    many = np.vstack((xyz[:3], np.arange(162).reshape(54, 3) / 100 + 0.2))
    term = HarmonicTarget(CoordinationDifferenceCV(0, 1, tuple(range(3, 57))), 0.0, 1.0)
    assert np.isfinite(evaluate_bias(many, [term]).gradient).all()


def test_callback_units_freezing_and_snapshot(xyz: np.ndarray, terms: tuple) -> None:
    def physical(positions: np.ndarray) -> tuple[float, np.ndarray]:
        return float((positions**2).sum() / HARTREE_TO_KCAL_MOL), 2 * positions

    mutable_h = [3, 4, 5]
    mutable_terms = [HarmonicTarget(CoordinationCV(0, mutable_h), 1.0, 3.0), *terms]
    callback = ReactionBiasCallback(physical, xyz, [2, 0, 4], mutable_terms)
    mutable_h.clear()
    mutable_terms.clear()
    xyz_saved = xyz.copy()
    xyz[:] += 1.0
    active = xyz_saved[[2, 0, 4]].ravel() + 0.2
    energy, gradient = callback(active)
    full = callback.coordinates(active)
    np.testing.assert_array_equal(full[[1, 3, 5]], xyz_saved[[1, 3, 5]])
    assert callback.last_evaluation.physical_energy == pytest.approx((full**2).sum())
    assert energy == pytest.approx(
        (full**2).sum() + callback.last_evaluation.bias.energy
    )
    np.testing.assert_allclose(
        gradient, callback.last_evaluation.gradient[[2, 0, 4]].ravel()
    )
    for bad in [np.zeros(4), np.full(9, np.nan), np.zeros((3, 3))]:
        with pytest.raises(ValueError):
            callback(bad)
        assert callback.last_evaluation is None


@pytest.mark.parametrize(
    "cv",
    [
        DistanceCV(0, 0),
        DistanceCV(-1, 2),
        DistanceCV(0, 8),
        DistanceCV(0.0, 2),
        CoordinationCV(0, ()),
        CoordinationCV(0, (3, 3)),
        CoordinationCV(0, (0, 3)),
        CoordinationCV(0, (3, 8)),
        CoordinationCV(0, (3.0, 4.0)),
        CoordinationCV(0, (3, 4), 0),
        CoordinationCV(0, (3, 4), np.nan),
        CoordinationDifferenceCV(0, 0, (3, 4)),
        object(),
    ],
)
def test_invalid_mapping(xyz: np.ndarray, cv: object) -> None:
    with pytest.raises(ValueError):
        evaluate_bias(xyz, [HarmonicTarget(cv, 1.0, 2.0)])


@pytest.mark.parametrize(
    "changes",
    [
        {"target": np.nan},
        {"spring": np.inf},
        {"spring": -1},
        {"target": 0},
        {"target": [1.0]},
        {"spring": [1.0]},
    ],
)
def test_invalid_target(xyz: np.ndarray, terms: tuple, changes: dict) -> None:
    with pytest.raises(ValueError):
        evaluate_bias(xyz, [replace(terms[0], **changes)])


@pytest.mark.parametrize("active", [[], [1, 1], [-1], [8], [1.5], [[1]], [True]])
def test_invalid_active(xyz: np.ndarray, terms: tuple, active: list) -> None:
    with pytest.raises(ValueError):
        ReactionBiasCallback(lambda x: (0.0, x), xyz, active, terms)


def test_other_invalid_inputs(xyz: np.ndarray, terms: tuple) -> None:
    for positions, targets in [
        (xyz, []),
        (xyz, [0]),
        (xyz[:0], terms),
        (xyz.ravel(), terms),
        (xyz * np.nan, terms),
    ]:
        with pytest.raises(ValueError):
            evaluate_bias(positions, targets)
    zero = xyz.copy()
    zero[1] = zero[0]
    with pytest.raises(ValueError):
        evaluate_bias(zero, terms)
    for energy, gradient in [
        (np.nan, xyz),
        ([0.0], xyz),
        (0.0, xyz[:2]),
        (0.0, xyz * np.nan),
    ]:
        with pytest.raises(ValueError):
            evaluate_biased(xyz, lambda x: (energy, gradient), terms)


def test_numpy_metadata(xyz: np.ndarray) -> None:
    term = HarmonicTarget(
        CoordinationCV(np.int64(0), np.array([3, 4]), np.float64(1.25)),
        np.float64(1.0),
        np.float64(2.0),
    )
    evaluate_bias(xyz, [term])
    assert json.loads(json.dumps(term.to_metadata()))["cv"]["center"] == 0


def test_linear_combination_gradient_units_and_serialization(xyz: np.ndarray) -> None:
    from nyl12_covalent.reaction_bias import LinearCombinationCV, WeightedCV

    hydrogens = [3, 4, 5]
    weighted = [
        WeightedCV(DistanceCV(0, 1), 1.0),
        WeightedCV(DistanceCV(1, 2), -1.0),
        WeightedCV(CoordinationDifferenceCV(2, 0, hydrogens), 1.2),
    ]
    cv = LinearCombinationCV(weighted, output_unit="angstrom")
    hydrogens.clear()
    weighted.clear()
    target = HarmonicTarget(cv, -0.3, 4.0)  # kcal/mol/Å²
    result = evaluate_bias(xyz, [target])
    parts = [
        evaluate_bias(xyz, [HarmonicTarget(term.cv, 1.0, 0.0)]).values[0]
        for term in cv.terms
    ]
    scalar = parts[0] - parts[1] + 1.2 * parts[2]
    assert result.values[0] == pytest.approx(scalar)
    assert result.energy == pytest.approx(0.5 * 4.0 * (scalar + 0.3) ** 2)
    numerical = np.empty_like(xyz)
    for index in np.ndindex(xyz.shape):
        plus, minus = xyz.copy(), xyz.copy()
        plus[index] += 1e-6
        minus[index] -= 1e-6
        numerical[index] = (
            evaluate_bias(plus, [target]).energy - evaluate_bias(minus, [target]).energy
        ) / 2e-6
    np.testing.assert_allclose(result.gradient, numerical, atol=3e-8)
    np.testing.assert_allclose(result.gradient.sum(axis=0), 0.0, atol=1e-14)
    metadata = json.loads(json.dumps(target.to_metadata()))
    assert metadata["cv"]["output_unit"] == "angstrom"
    assert metadata["cv"]["terms"][2]["weight"] == 1.2
    assert metadata["cv"]["terms"][2]["cv"]["hydrogens"] == [3, 4, 5]
    callback = ReactionBiasCallback(lambda x: (0.0, x * 0), xyz, [0, 1], [target])
    assert callback(xyz[[0, 1]].ravel())[0] == pytest.approx(result.energy)


def test_linear_combination_invalid() -> None:
    from nyl12_covalent.reaction_bias import LinearCombinationCV, WeightedCV

    for weight in [np.nan, np.inf, [1.0]]:
        with pytest.raises(ValueError):
            WeightedCV(DistanceCV(0, 1), weight)
    for terms, unit in [
        ((), "angstrom"),
        ((DistanceCV(0, 1),), "angstrom"),
        ((WeightedCV(DistanceCV(0, 1), 1.0),), ""),
        ((WeightedCV(DistanceCV(0, 1), 1.0),), None),
    ]:
        with pytest.raises(ValueError):
            LinearCombinationCV(terms, unit)
    nested = LinearCombinationCV((WeightedCV(DistanceCV(0, 1), 1.0),), "angstrom")
    with pytest.raises(ValueError):
        WeightedCV(nested, 1.0)
