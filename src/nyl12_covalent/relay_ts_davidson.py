"""Bounded all-active block-Davidson diagnostics with native Hessian products."""

from __future__ import annotations

from collections.abc import Callable
import logging
from pathlib import Path

import numpy as np

from .stationary_native import NativeState, digest, write_json

LOGGER = logging.getLogger(__name__)
Array = np.ndarray


def orthogonalize(vectors: Array, basis: Array) -> Array:
    """Return independent orthonormal columns after twice removing a basis."""
    vectors = vectors.copy()
    for _ in range(2):
        vectors -= basis @ (basis.T @ vectors)
    left, singular, _ = np.linalg.svd(vectors, full_matrices=False)
    if not singular.size:
        return vectors[:, :0]
    return left[:, singular > 1e-10]


def davidson(
    operator: Callable[[Array], Array],
    precondition: Callable[[Array, float], Array],
    initial: Array,
    *,
    count: int = 3,
    budget: int = 74,
    tolerance: float = 1e-6,
) -> dict[str, object]:
    """Find low Ritz candidates within a strict operator-product budget.

    Args:
        operator: Full-space Hessian-vector product.
        precondition: Approximate inverse shifted Hessian for residual correction.
        initial: Full-space initial vectors including the chemical mode.
        count: Number of lowest Ritz candidates requested.
        budget: Strict maximum operator products, including initial vectors.
        tolerance: Full-space Ritz residual convergence tolerance.

    Returns:
        Low Ritz values/vectors, actual residuals and budget status. Finite
        iteration Ritz values alone do not certify complete Hessian inertia.
    """
    if (
        initial.ndim != 2
        or not np.isfinite(initial).all()
        or not 1 <= count < budget
        or tolerance <= 0
    ):
        raise ValueError("Invalid block eigensolver inputs")
    basis = orthogonalize(initial, np.empty((initial.shape[0], 0)))
    if not count <= basis.shape[1] <= budget:
        raise ValueError("Insufficient independent initial vectors or excessive block")
    products = np.column_stack([operator(v) for v in basis.T])
    while True:
        projected = basis.T @ products
        values, coefficients = np.linalg.eigh((projected + projected.T) / 2)
        values, coefficients = values[:count], coefficients[:, :count]
        modes = basis @ coefficients
        residual = products @ coefficients - modes * values[None, :]
        norms = np.linalg.norm(residual, axis=0)
        LOGGER.info(
            "Davidson %s products; low=%s residual=%s", basis.shape[1], values, norms
        )
        if np.max(norms) <= tolerance or basis.shape[1] >= budget:
            break
        corrections = np.column_stack(
            [
                precondition(-residual[:, i], float(values[i]))
                for i in range(count)
                if norms[i] > tolerance
            ]
        )
        additions = orthogonalize(corrections, basis)
        if not additions.shape[1]:
            additions = orthogonalize(residual, basis)
        additions = additions[:, : budget - basis.shape[1]]
        if not additions.shape[1]:
            break
        products = np.column_stack((products, *[operator(v) for v in additions.T]))
        basis = np.column_stack((basis, additions))
    return {
        "values": values,
        "modes": modes,
        "residuals": norms,
        "products": basis.shape[1],
        "residual_converged": bool(np.max(norms) <= tolerance),
        "projected_antisymmetry_max_abs": float(
            np.max(np.abs(projected - projected.T))
        ),
    }


class NativeHessianProducts:
    """Full-original-active-space central analytic-gradient Hessian products."""

    def __init__(
        self, state: NativeState, output: Path, step: float = 0.000125
    ) -> None:
        self.state = state
        self.output = output
        self.step = step
        self.count = 0

    def product(self, vector: Array, *, half: bool = False) -> Array:
        """Return a mass-weighted product with bounded maximum Cartesian step."""
        if vector.shape != self.state.inv_sqrt.shape or not np.isfinite(vector).all():
            raise ValueError("Invalid all-active mode inventory")
        cartesian = vector * self.state.inv_sqrt
        maximum = float(np.max(np.abs(cartesian)))
        if maximum <= 1e-14 or self.count >= 80:
            raise ValueError("Degenerate vector or exhausted80-product budget")
        delta = self.step / maximum / (2 if half else 1)
        directory = self.output / f"product_{self.count:03d}"
        directory.mkdir()
        np.save(directory / "mass_weighted_vector.npy", vector, allow_pickle=False)
        gradients = []
        for sign, label in ((1, "plus"), (-1, "minus")):
            xyz = self.state.xyz.copy()
            xyz[self.state.active] += sign * delta * cartesian.reshape((-1, 3))
            _, gradient = self.state.gradient(xyz, directory / label)
            gradients.append(gradient.reshape(-1))
        result = (gradients[0] - gradients[1]) / (2 * delta) * self.state.inv_sqrt
        np.save(directory / "product.npy", result, allow_pickle=False)
        write_json(
            directory / "specification.json",
            {
                "scalar_displacement": delta,
                "maximum_Cartesian_displacement_A": self.step / (2 if half else 1),
                "half_step": half,
                "active_atoms": 4165,
                "units": "Eh/(A2 amu)",
                "vector_sha256": digest(directory / "mass_weighted_vector.npy"),
                "product_sha256": digest(directory / "product.npy"),
            },
        )
        self.count += 1
        write_json(
            self.output / "progress.json",
            {"operator_products": self.count, "maximum": 80},
        )
        return result


def run_davidson(
    model: Path, restart: Path, screen: Path, output: Path
) -> dict[str, object]:
    """Run at most80all-active products, retaining independent half-step probes.

    Args:
        model: Unchanged original model and native precision library.
        restart: Original force-accepted TS candidate.
        screen: Existing partial-QM derivative screen and chemical mode.
        output: Fresh all-active diagnostic directory.

    Returns:
        Full-active Ritz/residual diagnostics with unresolved signs explicit.
        A bounded subspace calculation is not complete full-Hessian inertia.
    """
    output.mkdir(parents=True, exist_ok=False)
    state = NativeState(model, restart)
    columns = np.load(screen / "partial_columns.npy")
    raw = np.load(screen / "hessian_raw.npy")
    partial = (raw + raw.T) / 2
    partial_values, partial_modes = np.linalg.eigh(partial)
    stiffness = float(np.median(np.diag(partial) / state.inv_sqrt[columns] ** 2))
    if stiffness <= 0:
        raise ValueError("Nonpositive partial-block median diagonal preconditioner")
    diagonal = stiffness * state.inv_sqrt**2

    def precondition(vector: Array, value: float) -> Array:
        result = vector / np.maximum(diagonal - value, 1e-5)
        denom = np.maximum(np.abs(partial_values - value), 1e-5)
        result[columns] = partial_modes @ ((partial_modes.T @ vector[columns]) / denom)
        return result

    chemical = np.load(screen / "lifted_cartesian_direction.npy") / state.inv_sqrt
    chemical /= np.linalg.norm(chemical)
    rng = np.random.default_rng(20261125)
    initial = np.column_stack((chemical, rng.normal(size=(len(chemical), 3))))
    write_json(
        output / "specification.json",
        {
            "restart_sha256": digest(restart),
            "partial_hessian_sha256": digest(screen / "hessian_raw.npy"),
            "active_atoms": 4165,
            "dimension": len(chemical),
            "maximum_HVP": 80,
            "initial_vectors": "Lifted chemical mode plus3deterministic shell-inclusive Gaussian vectors",
            "seed": 20261125,
            "maximum_Cartesian_step_A": 0.000125,
            "preconditioner": "Partial QM eigensystem; positive median QM stiffness times inverse masses for MM; algorithm only, unchanged physical Hessian",
            "median_unweighted_QM_diagonal": stiffness,
            "residual_tolerance": 1e-6,
            "full_active_index_verified": False,
        },
    )
    native = NativeHessianProducts(state, output)
    result = davidson(native.product, precondition, initial)
    values, modes = result["values"], result["modes"]
    np.save(output / "ritz_modes.npy", modes, allow_pickle=False)
    probes = []
    for i, mode in enumerate(modes.T):
        coarse = native.product(mode)
        fine = native.product(mode, half=True)
        rayleigh = [float(mode @ coarse), float(mode @ fine)]
        probes.append(
            {
                "mode": i,
                "Rayleigh_coarse_and_half": rayleigh,
                "product_half_step_difference_norm": float(
                    np.linalg.norm(coarse - fine)
                ),
                "independent_fine_residual_norm": float(
                    np.linalg.norm(fine - values[i] * mode)
                ),
            }
        )
    receipt = {
        "lowest_Ritz_values": values.tolist(),
        "full_space_Ritz_residual_norms": result["residuals"].tolist(),
        "residual_converged": result["residual_converged"],
        "operator_products": native.count,
        "projected_antisymmetry_max_abs": result["projected_antisymmetry_max_abs"],
        "independent_mode_probes": probes,
        "full_active_index_verified": False,
        "interpretation": "Bounded full-active low-mode diagnostic. Finite subspace Ritz/residual evidence and directional half-step probes do not establish complete inertia or an operator-norm error bound. Soft signs overlapping measured errors remain unresolved.",
    }
    write_json(output / "scientific_acceptance.json", receipt)
    return receipt
