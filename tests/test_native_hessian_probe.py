"""No-engine tests for finite-difference Hessian kernel and eigensolver."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse.linalg import ArpackNoConvergence

from nyl12_covalent import native_hessian_probe as probe
from nyl12_covalent.native_hessian_probe import (
    estimate_low_modes,
    finite_difference_hvp,
)


def test_hvp_matches_diagonal_mass_weighted_hessian() -> None:
    coordinates = np.zeros((3, 3))
    active = np.array([0, 2])
    masses = np.array([1.0, 4.0])
    inverse_sqrt = 1.0 / np.sqrt(masses)
    diagonal = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0])

    def gradient(x: np.ndarray) -> np.ndarray:
        full = np.zeros_like(x)
        full[0] = x[0] * np.array([2, 3, 4])
        full[2] = x[2] * np.array([5, 6, 7])
        return full

    vector = np.arange(1.0, 7.0)
    result = finite_difference_hvp(
        vector, coordinates, active, inverse_sqrt, gradient, 1e-4
    )
    expected = diagonal * inverse_sqrt.repeat(3) ** 2 * vector
    assert result == pytest.approx(expected)


def test_eigsh_finds_two_smallest_eigenvalues_and_residuals() -> None:
    matrix = np.diag([-0.4, 0.2, 0.5, 1.0, 2.0, 3.0])
    result = estimate_low_modes(
        lambda vector: matrix @ vector,
        dimension=6,
        count=2,
        tolerance=1e-10,
        max_iterations=100,
        subspace_size=5,
    )
    assert result["eigenvalues_hartree_per_angstrom2_amu"] == pytest.approx([-0.4, 0.2])
    assert max(result["residual_norms_hartree_per_angstrom2_amu"]) < 1e-8


def test_hvp_rejects_zero_direction() -> None:
    with pytest.raises(ValueError, match="nonzero"):
        finite_difference_hvp(
            np.zeros(3),
            np.zeros((1, 3)),
            np.array([0]),
            np.ones(1),
            lambda coords: coords,
            0.01,
        )


def test_eigsh_receives_normalized_physical_start_vector(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: dict[str, np.ndarray] = {}
    matrix = np.diag([-1.0, 0.5, 1.0, 2.0])

    def fake_eigsh(
        _operator: object, **kwargs: object
    ) -> tuple[np.ndarray, np.ndarray]:
        observed["v0"] = np.asarray(kwargs["v0"])
        vectors = np.eye(4)[:, :2]
        return np.array([-1.0, 0.5]), vectors

    monkeypatch.setattr(probe, "eigsh", fake_eigsh)
    initial = np.array([0.0, 0.0, 3.0, 4.0])
    result = estimate_low_modes(
        lambda vector: matrix @ vector,
        4,
        initial_vector=initial,
        subspace_size=4,
    )

    assert observed["v0"] == pytest.approx([0.0, 0.0, 0.6, 0.8])
    assert result["status"] == "converged"


def test_arpack_nonconvergence_retains_partial_eigenpairs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    partial_vectors = np.eye(4)[:, :1]
    error = ArpackNoConvergence("iteration cap", np.array([-0.25]), partial_vectors)

    def fail_eigsh(*_args: object, **_kwargs: object) -> None:
        raise error

    monkeypatch.setattr(probe, "eigsh", fail_eigsh)
    result = estimate_low_modes(
        lambda _vector: pytest.fail(
            "partial eigenpairs should not trigger residual HVP"
        ),
        4,
        initial_vector=np.ones(4),
        subspace_size=4,
    )

    assert result["status"] == "partial_convergence"
    assert result["eigenvalues_hartree_per_angstrom2_amu"] == [-0.25]
    assert result["eigenvectors"].shape == (4, 1)
    assert result["residual_norms_hartree_per_angstrom2_amu"] == []
