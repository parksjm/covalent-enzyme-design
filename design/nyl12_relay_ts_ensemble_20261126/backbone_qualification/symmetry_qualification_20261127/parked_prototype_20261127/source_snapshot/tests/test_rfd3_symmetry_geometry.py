"""Tests of independent D2 operators and coordinate residuals."""

import numpy as np
import pytest

from nyl12_covalent.rfd3_symmetry_geometry import (
    D2Operators,
    displacement_metrics,
    fit_d2_operators,
    measure_d2_residual,
    operator_closure_metrics,
    reference_projection,
    validate_correspondence,
)


def fixture_geometry() -> tuple[
    np.ndarray, list[list[tuple[str, int, str]]], D2Operators
]:
    """Create a translated, rotated exact D2 coupled alpha/beta assembly."""
    rng = np.random.default_rng(27)
    frame, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rotations = np.stack(
        [np.eye(3)] + [2 * np.outer(axis, axis) - np.eye(3) for axis in frame.T]
    )
    operators = D2Operators(rotations, np.array([11.0, -7.0, 19.0]))
    repeat = rng.normal(size=(12, 3)) + [5, 7, 9]
    xyz = reference_projection(np.stack([repeat] * 4), operators)
    keys = [("alpha" if i < 8 else "beta", i + 1, "CA") for i in range(12)]
    return xyz, [keys.copy() for _ in range(4)], operators


def test_fit_recovers_shared_operators_and_closure() -> None:
    xyz, keys, expected = fixture_geometry()
    actual = fit_d2_operators(xyz, keys)
    np.testing.assert_allclose(actual.rotations, expected.rotations, atol=1e-12)
    np.testing.assert_allclose(actual.center, expected.center, atol=1e-12)
    assert max(operator_closure_metrics(actual).values()) < 1e-12
    metrics = measure_d2_residual(
        xyz, actual, rmsd_limit_A=1e-8, max_deviation_limit_A=1e-8
    )
    assert metrics["within_tolerances"]
    assert all(np.linalg.det(rot) > 0.999999 for rot in actual.rotations)
    assert not actual.rotations.flags.writeable


def test_asymmetric_beta_translation_is_not_hidden_by_refitting() -> None:
    xyz, _, operators = fixture_geometry()
    xyz[2, 8:, 0] += 2
    metrics = measure_d2_residual(
        xyz, operators, rmsd_limit_A=0.01, max_deviation_limit_A=0.01
    )
    assert not metrics["within_tolerances"]
    assert metrics["backbone_max_deviation_A"] == pytest.approx(1.5)
    assert metrics["per_repeat"][2]["max_deviation_A"] == pytest.approx(1.5)


def test_projection_retains_reference_and_reports_displacement() -> None:
    xyz, _, operators = fixture_geometry()
    xyz[1, 3] += 1
    projected = reference_projection(xyz, operators)
    np.testing.assert_allclose(projected[0], xyz[0])
    assert displacement_metrics(xyz, projected)["max_deviation_A"] == pytest.approx(
        np.sqrt(3)
    )
    fixed = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    assert displacement_metrics(fixed, fixed)["rmsd_A"] == 0
    assert displacement_metrics(fixed, fixed + [0, 0, 0.5])["max_deviation_A"] == 0.5


@pytest.mark.parametrize(
    "poison", ["order", "duplicate", "segment", "empty", "residue", "missing"]
)
def test_poisoned_correspondence_rejected(poison: str) -> None:
    _, keys, _ = fixture_geometry()
    if poison == "order":
        keys[1].reverse()
    elif poison == "duplicate":
        keys[0][1] = keys[0][0]
    elif poison == "segment":
        keys = [[("alpha", key[1], key[2]) for key in repeat] for repeat in keys]
    elif poison == "empty":
        keys[0][0] = ("alpha", 1, "")
    elif poison == "residue":
        keys[0][0] = ("alpha", 0, "CA")
    else:
        keys.pop()
    with pytest.raises(ValueError):
        validate_correspondence(keys, 12)


@pytest.mark.parametrize("poison", ["nan", "shape", "rank", "axes"])
def test_invalid_fit_rejected(poison: str) -> None:
    xyz, keys, _ = fixture_geometry()
    if poison == "nan":
        xyz[0, 0, 0] = np.nan
    elif poison == "shape":
        xyz = xyz[:3]
    elif poison == "rank":
        xyz[:, :, 2] = 0
    else:
        xyz = np.stack([xyz[0]] * 4)
    with pytest.raises(ValueError):
        fit_d2_operators(xyz, keys)


@pytest.mark.parametrize(
    "poison",
    ["identity", "reflection", "nonorthogonal", "duplicate", "center", "quarterturn"],
)
def test_invalid_operators_rejected(poison: str) -> None:
    _, _, operators = fixture_geometry()
    rotations, center = operators.rotations.copy(), operators.center.copy()
    if poison == "identity":
        rotations[0] = rotations[1]
    elif poison == "reflection":
        rotations[1] *= -1
    elif poison == "nonorthogonal":
        rotations[1, 0, 0] += 0.2
    elif poison == "duplicate":
        rotations[2] = rotations[1]
    elif poison == "quarterturn":
        rotations[1] = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    else:
        center[0] = np.nan
    with pytest.raises(ValueError):
        D2Operators(rotations, center)


@pytest.mark.parametrize("limit", [0.0, -1.0, np.nan, np.inf])
def test_invalid_tolerance_rejected(limit: float) -> None:
    xyz, _, operators = fixture_geometry()
    with pytest.raises(ValueError):
        measure_d2_residual(
            xyz, operators, rmsd_limit_A=limit, max_deviation_limit_A=1.0
        )


@pytest.mark.parametrize(
    "source, output",
    [
        ([], []),
        ([[1.0, 2.0, 3.0]], [[np.nan, 2.0, 3.0]]),
        ([[1.0, 2.0, 3.0]], [[1.0, 2.0]]),
        ([[np.inf, 2.0, 3.0]], [[1.0, 2.0, 3.0]]),
    ],
)
def test_invalid_fixed_context_rejected(source: list, output: list) -> None:
    with pytest.raises(ValueError):
        displacement_metrics(np.array(source), np.array(output))


def test_tied_displacements_preserve_asymmetric_source() -> None:
    from nyl12_covalent.rfd3_symmetry_geometry import measure_tied_displacement_residual

    source, _, operators = fixture_geometry()
    source[2, 3] += [1.0, 3.0, -2.0]
    unchanged = measure_tied_displacement_residual(
        source, source, operators, rmsd_limit_A=1e-8, max_deviation_limit_A=1e-8
    )
    assert unchanged["within_tolerances"]
    assert unchanged["tied_displacement_rmsd_A"] == 0
    absolute = measure_d2_residual(
        source, operators, rmsd_limit_A=1e-8, max_deviation_limit_A=1e-8
    )
    assert not absolute["within_tolerances"]
    common_delta = np.arange(36).reshape(12, 3) * 0.001
    tied_delta = np.stack([common_delta @ rot.T for rot in operators.rotations])
    output = source + tied_delta
    result = measure_tied_displacement_residual(
        source, output, operators, rmsd_limit_A=1e-8, max_deviation_limit_A=1e-8
    )
    assert result["within_tolerances"]
    output[1, 5, 1] += 0.1
    result = measure_tied_displacement_residual(
        source, output, operators, rmsd_limit_A=1e-8, max_deviation_limit_A=1e-8
    )
    assert not result["within_tolerances"]
    assert result["tied_displacement_max_deviation_A"] == pytest.approx(0.075)
    fixed_context = np.array([[1.0, 2.0, 3.0], [3.0, 4.0, 5.0]])
    assert displacement_metrics(fixed_context, fixed_context)["max_deviation_A"] == 0


def test_tied_displacement_invalid_output_rejected() -> None:
    from nyl12_covalent.rfd3_symmetry_geometry import measure_tied_displacement_residual

    source, _, operators = fixture_geometry()
    with pytest.raises(ValueError):
        measure_tied_displacement_residual(
            source, source[:3], operators, rmsd_limit_A=0.1, max_deviation_limit_A=0.1
        )
