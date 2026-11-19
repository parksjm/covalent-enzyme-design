"""Bounded block eigensolver for a matrix-free mass-weighted Hessian."""

from __future__ import annotations

import warnings
from typing import Callable

import numpy as np
from scipy.sparse.linalg import LinearOperator, lobpcg


def estimate_block_low_modes(
    matvec: Callable[[np.ndarray], np.ndarray],
    dimension: int,
    initial_vector: np.ndarray,
    *,
    initial_vectors: np.ndarray | None = None,
    count: int = 2,
    tolerance: float = 1.0e-4,
    max_iterations: int = 30,
    seed: int = 20261114,
) -> dict[str, object]:
    """Estimate low Hessian modes with a dimer-seeded block iteration.

    Args:
        matvec: Symmetric Hessian-vector product in mass-weighted coordinates.
        dimension: Number of active Cartesian coordinates.
        initial_vector: Saved chemical dimer mode in that coordinate order.
        initial_vectors: Optional matrix of saved Ritz vectors, shape
            (dimension, at most count).
        count: Number of lowest algebraic modes sought.
        tolerance: Absolute eigenpair residual target.
        max_iterations: Hard limit on block iterations.
        seed: Reproducible seed for other initial block vectors.

    Returns:
        Eigenvalues, vectors, residuals and convergence diagnostics.

    Raises:
        ValueError: Dimensions, vector or solver controls are invalid.
    """
    if dimension < 5 * count or count < 1:
        raise ValueError("Block eigenproblem is too small for requested modes")
    if not 0.0 < tolerance < 1.0 or max_iterations < 1:
        raise ValueError("Invalid block solver controls")
    first = np.asarray(initial_vector, dtype=float).reshape(-1)
    if first.shape != (dimension,) or not np.isfinite(first).all():
        raise ValueError("Initial mode has invalid dimension or values")
    norm = np.linalg.norm(first)
    if norm <= 1.0e-12:
        raise ValueError("Initial dimer mode is zero")
    rng = np.random.default_rng(seed)
    if initial_vectors is None:
        block = np.column_stack([first / norm, rng.normal(size=(dimension, count - 1))])
    else:
        warm = np.asarray(initial_vectors, dtype=float)
        if (
            warm.ndim != 2
            or warm.shape[0] != dimension
            or not 1 <= warm.shape[1] <= count
        ):
            raise ValueError("Warm-start modes must match (dimension, at most count)")
        if not np.isfinite(warm).all() or np.any(
            np.linalg.norm(warm, axis=0) <= 1.0e-12
        ):
            raise ValueError("Warm-start modes must be finite and nonzero")
        if np.linalg.matrix_rank(warm) != warm.shape[1]:
            raise ValueError("Warm-start modes must be linearly independent")
        block = (
            warm.copy()
            if warm.shape[1] == count
            else np.column_stack([first / norm, warm])
        )
        if block.shape[1] < count:
            block = np.column_stack(
                [block, rng.normal(size=(dimension, count - block.shape[1]))]
            )
    block, _ = np.linalg.qr(block)

    def apply(vector: np.ndarray) -> np.ndarray:
        return np.asarray(matvec(np.asarray(vector).reshape(-1))).reshape(-1)

    operator = LinearOperator((dimension, dimension), matvec=apply, dtype=float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        values, vectors, history = lobpcg(
            operator,
            block,
            largest=False,
            tol=tolerance,
            maxiter=max_iterations,
            retResidualNormsHistory=True,
        )
    order = np.argsort(values)
    values, vectors = values[order], vectors[:, order]
    residuals = [
        float(np.linalg.norm(apply(vectors[:, i]) - values[i] * vectors[:, i]))
        for i in range(count)
    ]
    return {
        "status": (
            "converged"
            if all(residual <= tolerance for residual in residuals)
            else "partial_convergence"
        ),
        "eigenvalues_hartree_per_angstrom2_amu": values.tolist(),
        "eigenvectors": vectors,
        "residual_norms_hartree_per_angstrom2_amu": residuals,
        "iteration_residual_history": [
            np.asarray(entry, dtype=float).tolist() for entry in history
        ],
        "warnings": [str(item.message) for item in caught],
    }


def cache_hvp_direction(
    basis: list[np.ndarray],
    products: list[np.ndarray],
    vector: np.ndarray,
    product: np.ndarray,
) -> None:
    """Add an independent HVP direction and its image to an orthonormal cache.

    Args:
        basis: Mutable orthonormal basis vectors in mass-weighted coordinates.
        products: Mutable Hessian images corresponding to ``basis``.
        vector: Input vector for a completed Hessian-vector product.
        product: Hessian-vector product corresponding to ``vector``.

    Raises:
        ValueError: Vector or product is invalid or has the wrong dimension.
    """
    candidate = np.asarray(vector, dtype=float).reshape(-1).copy()
    image = np.asarray(product, dtype=float).reshape(-1).copy()
    if candidate.shape != image.shape or not np.isfinite(candidate).all():
        raise ValueError("Cached HVP vector and product must be finite and aligned")
    if not np.isfinite(image).all():
        raise ValueError("Cached HVP product must be finite")
    source_norm = float(np.linalg.norm(candidate))
    if source_norm <= 1.0e-12:
        return
    if basis:
        q = np.column_stack(basis)
        hq = np.column_stack(products)
        for _ in range(2):
            coefficients = q.T @ candidate
            candidate -= q @ coefficients
            image -= hq @ coefficients
    norm = float(np.linalg.norm(candidate))
    if norm > 1.0e-10 * source_norm:
        basis.append(candidate / norm)
        products.append(image / norm)


def project_cached_hvp_subspace(
    basis: list[np.ndarray],
    products: list[np.ndarray],
    *,
    count: int,
    tolerance: float,
) -> dict[str, object]:
    """Return lowest Ritz pairs and true residuals from cached HVPs.

    The cache must contain an orthonormal basis and its corresponding Hessian
    images. This projection preserves useful solver progress if an external HVP
    budget interrupts an iterative eigensolver.

    Args:
        basis: Orthonormal directions with completed Hessian-vector products.
        products: Hessian images corresponding to each basis direction.
        count: Number of lowest algebraic Ritz pairs to return.
        tolerance: Absolute eigenpair residual acceptance threshold.

    Returns:
        Partial or converged Ritz pairs, their true residuals, and subspace size.

    Raises:
        ValueError: The cache is inconsistent or too small for ``count`` modes.
    """
    if count < 1 or len(basis) != len(products) or len(basis) < count:
        raise ValueError("HVP cache must contain at least count paired directions")
    q = np.column_stack(basis)
    hq = np.column_stack(products)
    if q.shape != hq.shape or not np.isfinite(q).all() or not np.isfinite(hq).all():
        raise ValueError("HVP cache arrays must be finite and dimension matched")
    projected = q.T @ hq
    projected = 0.5 * (projected + projected.T)
    values, rotation = np.linalg.eigh(projected)
    order = np.argsort(values)[:count]
    values = values[order]
    vectors = q @ rotation[:, order]
    images = hq @ rotation[:, order]
    residuals = np.linalg.norm(images - vectors * values[None, :], axis=0)
    return {
        "status": (
            "converged" if np.all(residuals <= tolerance) else "partial_convergence"
        ),
        "eigenvalues_hartree_per_angstrom2_amu": values.tolist(),
        "eigenvectors": vectors,
        "residual_norms_hartree_per_angstrom2_amu": residuals.tolist(),
        "cached_subspace_dimension": q.shape[1],
    }
