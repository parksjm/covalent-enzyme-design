"""Assess local curvature on a declared nonlinear-constraint manifold."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from numpy.typing import NDArray

LOGGER = logging.getLogger(__name__)
FloatArray = NDArray[np.float64]


def _rank(values: FloatArray, scale: float) -> int:
    """Determine numerical rank using a declared relative singular-value scale."""
    return int(np.count_nonzero(values > scale * values[0])) if values.size else 0


def reduced_lagrangian_curvature(
    hessian: FloatArray,
    gradient: FloatArray,
    constraint_jacobian: FloatArray,
    constraint_hessians: FloatArray,
    *,
    excluded_directions: FloatArray,
    stationarity_tolerance: float,
    curvature_tolerance: float,
) -> dict[str, Any]:
    """Test stationarity and reduced curvature for equality constraints.

    Use ``L = E - lambda . g`` with ``g = 0`` and solve
    ``J.T lambda = gradient``. Project the Hessian of L, including second
    derivatives of nonlinear constraints, onto the remaining tangent space.
    This is not a mass-weighted vibrational analysis or force-field validation.

    Args:
        hessian: Symmetric objective Hessian, shape (n, n), Hartree/bohr^2.
        gradient: Objective gradient, shape (n,), Hartree/bohr.
        constraint_jacobian: Constraint Jacobian, shape (m, n).
        constraint_hessians: Constraint second derivatives, shape (m, n, n).
        excluded_directions: Independent Cartesian directions to remove, such
            as rigid motions, shape (n, k). Must lie in the constraint tangent.
        stationarity_tolerance: Maximum absolute residual Lagrangian gradient.
        curvature_tolerance: Eigenvalue magnitude below which sign is unresolved.

    Returns:
        Explicit constrained-stability assessment and unweighted eigenvalues.

    Raises:
        ValueError: If shapes, symmetry, finite values, tolerances or ranks fail.
    """
    hessian = np.asarray(hessian, dtype=float)
    gradient = np.asarray(gradient, dtype=float)
    jacobian = np.asarray(constraint_jacobian, dtype=float)
    second = np.asarray(constraint_hessians, dtype=float)
    removed = np.asarray(excluded_directions, dtype=float)
    dimension = gradient.size
    if gradient.shape != (dimension,) or dimension < 2:
        raise ValueError("Expected a one-dimensional Cartesian gradient")
    if hessian.shape != (dimension, dimension):
        raise ValueError("Objective Hessian has the wrong shape")
    if (
        jacobian.ndim != 2
        or jacobian.shape[1] != dimension
        or not 0 < jacobian.shape[0] < dimension
    ):
        raise ValueError("Expected one or more independent constraint rows")
    constraints = jacobian.shape[0]
    if second.shape != (constraints, dimension, dimension):
        raise ValueError("Constraint Hessians have the wrong shape")
    if removed.ndim != 2 or removed.shape[0] != dimension:
        raise ValueError("Excluded directions have the wrong shape")
    for value in (hessian, gradient, jacobian, second, removed):
        if not np.isfinite(value).all():
            raise ValueError("All curvature inputs must be finite")
    for value in (stationarity_tolerance, curvature_tolerance):
        if not np.isfinite(value) or value <= 0:
            raise ValueError("Both assessment tolerances must be finite and positive")
    if not np.allclose(hessian, hessian.T, rtol=1e-10, atol=1e-12):
        raise ValueError("Objective Hessian is not symmetric")
    if not np.allclose(second, second.transpose(0, 2, 1), rtol=1e-10, atol=1e-12):
        raise ValueError("Constraint Hessians are not symmetric")
    rank_scale = 1e-10
    _, singular, _ = np.linalg.svd(jacobian, full_matrices=False)
    if _rank(singular, rank_scale) != constraints:
        raise ValueError("Constraint Jacobian has dependent rows")
    if removed.shape[1]:
        rigid_basis, singular, _ = np.linalg.svd(removed, full_matrices=False)
        if _rank(singular, rank_scale) != removed.shape[1]:
            raise ValueError("Excluded directions are dependent")
        rigid_basis = rigid_basis[:, : removed.shape[1]]
        if np.linalg.norm(jacobian @ rigid_basis) > 1e-6 * np.linalg.norm(jacobian):
            raise ValueError("Excluded directions do not lie in the constraint tangent")
    else:
        rigid_basis = np.empty((dimension, 0))
    # Row normalization makes the joint rank test independent of constraint units.
    rows = np.vstack(
        (jacobian / np.linalg.norm(jacobian, axis=1)[:, None], rigid_basis.T)
    )
    _, singular, right = np.linalg.svd(rows, full_matrices=True)
    rank = _rank(singular, rank_scale)
    if rank != constraints + removed.shape[1]:
        raise ValueError("Constraints and excluded directions are not independent")
    tangent = right[rank:].T
    if tangent.shape[1] == 0:
        raise ValueError("No internal tangent coordinates remain")
    multipliers = np.linalg.lstsq(jacobian.T, gradient, rcond=rank_scale)[0]
    residual = gradient - jacobian.T @ multipliers
    lagrangian_hessian = hessian - np.einsum("i,ijk->jk", multipliers, second)
    reduced = tangent.T @ lagrangian_hessian @ tangent
    eigenvalues = np.linalg.eigvalsh((reduced + reduced.T) / 2)
    minimum = float(eigenvalues[0])
    classification = (
        "negative"
        if minimum < -curvature_tolerance
        else "positive_definite"
        if minimum > curvature_tolerance
        else "unresolved"
    )
    residual_max = float(np.max(np.abs(residual)))
    stationary = residual_max <= stationarity_tolerance
    accepted = bool(stationary and classification == "positive_definite")
    LOGGER.info(
        "Constrained curvature %s; stationary=%s; accepted=%s",
        classification,
        stationary,
        accepted,
    )
    return {
        "accepted": accepted,
        "scope": "Local stability only under declared equality constraints",
        "lagrangian_convention": "L = E - lambda . g; J.T lambda = gradient",
        "lagrange_multipliers": multipliers.tolist(),
        "gradient_residual_max": residual_max,
        "stationarity_tolerance": stationarity_tolerance,
        "stationarity_passed": stationary,
        "curvature_classification": classification,
        "curvature_tolerance": curvature_tolerance,
        "reduced_eigenvalues": eigenvalues.tolist(),
        "eigenvalue_units": "Hartree/bohr^2; orthonormal Cartesian tangent; not mass weighted",
        "tangent_dimension": tangent.shape[1],
        "constraint_rank": constraints,
        "excluded_direction_count": removed.shape[1],
        "forcefield_validated": False,
        "unconstrained_minimum_claimed": False,
    }
