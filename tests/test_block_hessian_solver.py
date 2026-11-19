"""Checks for the dimer-seeded block Hessian solver."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.block_hessian_solver import (
    cache_hvp_direction,
    estimate_block_low_modes,
    project_cached_hvp_subspace,
)


def test_recovers_two_negative_modes_on_known_hessian() -> None:
    """A second negative mode must be visible rather than called first order."""
    diagonal = np.concatenate(([-4.0, -1.0], np.arange(1.0, 39.0)))
    initial = np.zeros(len(diagonal))
    initial[0] = 1.0
    result = estimate_block_low_modes(
        lambda vector: diagonal * vector,
        len(diagonal),
        initial,
        tolerance=1.0e-6,
        max_iterations=80,
    )
    assert result["status"] == "converged"
    np.testing.assert_allclose(
        result["eigenvalues_hartree_per_angstrom2_amu"],
        [-4.0, -1.0],
        atol=1.0e-6,
    )
    assert max(result["residual_norms_hartree_per_angstrom2_amu"]) < 1.0e-6


def test_rejects_zero_initial_mode() -> None:
    """An absent dimer mode cannot silently become a random TS search."""
    with pytest.raises(ValueError, match="zero"):
        estimate_block_low_modes(
            lambda vector: vector,
            20,
            np.zeros(20),
        )


def test_reuses_saved_ritz_subspace() -> None:
    """Saved Ritz vectors seed the next block without random replacement."""
    diagonal = np.concatenate(([-4.0, -1.0], np.arange(1.0, 39.0)))
    warm = np.eye(len(diagonal))[:, :2]
    result = estimate_block_low_modes(
        lambda vector: diagonal * vector,
        len(diagonal),
        np.ones(len(diagonal)),
        initial_vectors=warm,
        tolerance=1.0e-6,
        max_iterations=10,
    )

    assert result["status"] == "converged"
    np.testing.assert_allclose(
        result["eigenvalues_hartree_per_angstrom2_amu"],
        [-4.0, -1.0],
        atol=1.0e-6,
    )


def test_rejects_mismatched_ritz_subspace() -> None:
    """Warm-start shape errors fail before any Hessian products."""
    with pytest.raises(ValueError, match="match"):
        estimate_block_low_modes(
            lambda vector: vector,
            20,
            np.ones(20),
            initial_vectors=np.ones((20, 3)),
        )


def test_wider_block_resolves_clustered_low_modes() -> None:
    """A wider block resolves a cluster missed by a two-vector start."""
    low = np.array([-0.4, -0.3, -0.2, 0.01, 0.02, 0.03])
    diagonal = np.concatenate((low, np.linspace(0.1, 4.0, 54)))
    warm = np.eye(len(diagonal))[:, :2]
    result = estimate_block_low_modes(
        lambda vector: diagonal * vector,
        len(diagonal),
        np.eye(len(diagonal))[:, 2],
        count=6,
        initial_vectors=warm,
        tolerance=1.0e-6,
        max_iterations=80,
    )

    assert result["status"] == "converged"
    np.testing.assert_allclose(
        result["eigenvalues_hartree_per_angstrom2_amu"], low, atol=1.0e-6
    )


def test_hvp_cache_recovers_ritz_pairs_after_solver_budget_stop() -> None:
    """Cached products retain exact low modes when an outer cap interrupts."""
    diagonal = np.array([-3.0, -1.0, 2.0, 4.0])
    basis: list[np.ndarray] = []
    products: list[np.ndarray] = []
    for vector in np.eye(len(diagonal))[:, :3].T:
        cache_hvp_direction(basis, products, vector, diagonal * vector)

    result = project_cached_hvp_subspace(basis, products, count=2, tolerance=1.0e-10)

    assert result["status"] == "converged"
    assert result["cached_subspace_dimension"] == 3
    np.testing.assert_allclose(
        result["eigenvalues_hartree_per_angstrom2_amu"], [-3.0, -1.0]
    )
    np.testing.assert_allclose(
        result["residual_norms_hartree_per_angstrom2_amu"], [0.0, 0.0], atol=1e-12
    )


def test_cached_ritz_residual_uses_full_hvp_image() -> None:
    """Projection reports unresolved residual outside the cached subspace."""
    matrix = np.diag([-1.0, 1.0, 3.0])
    direction = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    basis: list[np.ndarray] = []
    products: list[np.ndarray] = []
    cache_hvp_direction(basis, products, direction, matrix @ direction)
    result = project_cached_hvp_subspace(basis, products, count=1, tolerance=1.0e-8)

    assert result["status"] == "partial_convergence"
    assert result["eigenvalues_hartree_per_angstrom2_amu"] == pytest.approx([0.0])
    assert result["residual_norms_hartree_per_angstrom2_amu"] == pytest.approx([1.0])
