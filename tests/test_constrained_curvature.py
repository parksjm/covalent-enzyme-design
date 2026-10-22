"""Analytic regression cases for constrained Lagrangian curvature."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.constrained_curvature import reduced_lagrangian_curvature


TOLERANCE = 1e-8


def _circle_case(point: tuple[float, float], gradient: np.ndarray) -> dict[str, object]:
    """Evaluate E=x on the unit circle at a supplied point."""
    x, y = point
    del y  # The gradient is supplied explicitly to test stationarity handling.
    return reduced_lagrangian_curvature(
        hessian=np.zeros((2, 2)),
        gradient=gradient,
        constraint_jacobian=np.array([[x, point[1]]]),
        constraint_hessians=np.array([np.eye(2)]),
        excluded_directions=np.empty((2, 0)),
        stationarity_tolerance=TOLERANCE,
        curvature_tolerance=TOLERANCE,
    )


def test_circle_minimum_has_positive_lagrangian_curvature() -> None:
    """E=x on x²+y²=1 has a constrained minimum at (-1, 0)."""
    result = _circle_case((-1.0, 0.0), np.array([1.0, 0.0]))

    assert result["accepted"] is True
    assert result["stationarity_passed"] is True
    assert result["curvature_classification"] == "positive_definite"
    assert result["lagrange_multipliers"] == pytest.approx([-1.0])
    assert result["reduced_eigenvalues"] == pytest.approx([1.0])
    assert result["tangent_dimension"] == 1
    assert result["forcefield_validated"] is False
    assert result["unconstrained_minimum_claimed"] is False


def test_circle_maximum_has_negative_lagrangian_curvature() -> None:
    """E=x on x²+y²=1 has a constrained maximum at (+1, 0)."""
    result = _circle_case((1.0, 0.0), np.array([1.0, 0.0]))

    assert result["accepted"] is False
    assert result["curvature_classification"] == "negative"
    assert result["lagrange_multipliers"] == pytest.approx([1.0])
    assert result["reduced_eigenvalues"] == pytest.approx([-1.0])


def test_nonstationary_circle_point_is_not_accepted() -> None:
    """A tangent gradient component must fail the stationarity check."""
    result = _circle_case((-1.0, 0.0), np.array([1.0, 0.1]))

    assert result["accepted"] is False
    assert result["stationarity_passed"] is False
    assert result["gradient_residual_max"] == pytest.approx(0.1)


def test_curvature_is_restricted_to_constraint_tangent() -> None:
    """A negative fixed-coordinate mode is absent from the feasible tangent."""
    result = reduced_lagrangian_curvature(
        hessian=np.diag([-100.0, 2.0, 3.0]),
        gradient=np.array([2.0, 0.0, 0.0]),
        constraint_jacobian=np.array([[1.0, 0.0, 0.0]]),
        constraint_hessians=np.zeros((1, 3, 3)),
        excluded_directions=np.empty((3, 0)),
        stationarity_tolerance=TOLERANCE,
        curvature_tolerance=TOLERANCE,
    )

    assert result["accepted"] is True
    assert result["reduced_eigenvalues"] == pytest.approx([2.0, 3.0])
    assert result["tangent_dimension"] == 2


def test_excluded_direction_is_removed_from_feasible_tangent() -> None:
    """An excluded z direction leaves only the y curvature eigenvalue."""
    result = reduced_lagrangian_curvature(
        hessian=np.diag([-100.0, 2.0, 3.0]),
        gradient=np.array([2.0, 0.0, 0.0]),
        constraint_jacobian=np.array([[1.0, 0.0, 0.0]]),
        constraint_hessians=np.zeros((1, 3, 3)),
        excluded_directions=np.array([[0.0], [0.0], [1.0]]),
        stationarity_tolerance=TOLERANCE,
        curvature_tolerance=TOLERANCE,
    )

    assert result["accepted"] is True
    assert result["reduced_eigenvalues"] == pytest.approx([2.0])
    assert result["tangent_dimension"] == 1


def test_near_zero_tangent_eigenvalue_is_unresolved() -> None:
    """A flat feasible mode cannot establish a strict constrained minimum."""
    result = reduced_lagrangian_curvature(
        hessian=np.diag([-100.0, 0.0, 2.0]),
        gradient=np.array([2.0, 0.0, 0.0]),
        constraint_jacobian=np.array([[1.0, 0.0, 0.0]]),
        constraint_hessians=np.zeros((1, 3, 3)),
        excluded_directions=np.empty((3, 0)),
        stationarity_tolerance=TOLERANCE,
        curvature_tolerance=TOLERANCE,
    )

    assert result["accepted"] is False
    assert result["curvature_classification"] == "unresolved"
    assert result["reduced_eigenvalues"] == pytest.approx([0.0, 2.0])


@pytest.mark.parametrize(
    ("hessian", "gradient", "jacobian", "constraint_hessians", "excluded"),
    [
        (
            np.array([[np.nan]]),
            np.array([0.0]),
            np.array([[1.0]]),
            np.zeros((1, 1, 1)),
            np.empty((1, 0)),
        ),
        (
            np.array([[1.0, 1.0], [0.0, 1.0]]),
            np.zeros(2),
            np.array([[1.0, 0.0]]),
            np.zeros((1, 2, 2)),
            np.empty((2, 0)),
        ),
        (
            np.eye(2),
            np.zeros(3),
            np.array([[1.0, 0.0]]),
            np.zeros((1, 2, 2)),
            np.empty((2, 0)),
        ),
    ],
)
def test_invalid_matrix_values_or_shapes_raise_value_error(
    hessian: np.ndarray,
    gradient: np.ndarray,
    jacobian: np.ndarray,
    constraint_hessians: np.ndarray,
    excluded: np.ndarray,
) -> None:
    """Malformed or non-finite curvature inputs fail explicitly."""
    with pytest.raises(ValueError):
        reduced_lagrangian_curvature(
            hessian=hessian,
            gradient=gradient,
            constraint_jacobian=jacobian,
            constraint_hessians=constraint_hessians,
            excluded_directions=excluded,
            stationarity_tolerance=TOLERANCE,
            curvature_tolerance=TOLERANCE,
        )


@pytest.mark.parametrize(
    ("jacobian", "excluded"),
    [
        (np.array([[1.0, 0.0], [2.0, 0.0]]), np.empty((2, 0))),
        (np.array([[1.0, 0.0]]), np.array([[0.0, 1.0], [0.0, 1.0]])),
        (np.array([[1.0, 0.0]]), np.array([[1.0], [0.0]])),
    ],
)
def test_dependent_or_incompatible_subspaces_raise_value_error(
    jacobian: np.ndarray, excluded: np.ndarray
) -> None:
    """Constraint rows and excluded directions must define valid subspaces."""
    n = 2
    with pytest.raises(ValueError):
        reduced_lagrangian_curvature(
            hessian=np.eye(n),
            gradient=np.zeros(n),
            constraint_jacobian=jacobian,
            constraint_hessians=np.zeros((jacobian.shape[0], n, n)),
            excluded_directions=excluded,
            stationarity_tolerance=TOLERANCE,
            curvature_tolerance=TOLERANCE,
        )


@pytest.mark.parametrize("tolerances", [(0.0, 1e-8), (1e-8, -1e-8)])
def test_nonpositive_tolerances_raise_value_error(
    tolerances: tuple[float, float],
) -> None:
    """Both numerical tolerances must be positive."""
    with pytest.raises(ValueError):
        reduced_lagrangian_curvature(
            hessian=np.eye(1),
            gradient=np.zeros(1),
            constraint_jacobian=np.array([[1.0]]),
            constraint_hessians=np.zeros((1, 1, 1)),
            excluded_directions=np.empty((1, 0)),
            stationarity_tolerance=tolerances[0],
            curvature_tolerance=tolerances[1],
        )
