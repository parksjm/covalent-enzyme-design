"""Numerical mass-weighted steepest-descent path steps.

This module provides a predictor-corrector pilot, not a formal IRC. A formal
IRC requires an initial Hessian eigenmode and a specified step/integration
protocol.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]
EnergyGradient = Callable[[FloatArray], tuple[float, FloatArray]]


@dataclass(frozen=True)
class DescentStep:
    """Result of one attempted mass-weighted predictor-corrector step."""

    coordinates: FloatArray
    energy: float
    step_size: float
    accepted: bool


def _finite_array(value: NDArray, name: str) -> None:
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must contain only finite values")


def _direction(
    gradient_x: FloatArray, masses: FloatArray, fixed_atoms: BoolArray
) -> FloatArray:
    """Return a unit downhill direction in mass-weighted coordinates."""
    gradient_q = gradient_x / np.sqrt(masses)[:, None]
    gradient_q[fixed_atoms] = 0.0
    norm = float(np.linalg.norm(gradient_q))
    if norm == 0.0:
        raise ValueError("Cannot determine descent direction from zero gradient")
    return -gradient_q / norm


def mass_weighted_descent_step(
    coordinates: FloatArray,
    masses: FloatArray,
    energy_gradient: EnergyGradient,
    step_size: float,
    fixed_atoms: BoolArray | None = None,
    initial_direction_q: FloatArray | None = None,
    minimum_step_size: float = 1e-6,
) -> DescentStep:
    """Attempt one mass-weighted Heun descent step with energy rejection.

    Args:
        coordinates: Cartesian coordinates with shape ``(n_atoms, 3)``.
        masses: Positive atomic masses with shape ``(n_atoms,)``.
        energy_gradient: Pure evaluator returning ``(energy, dE/dx)``.
        step_size: Positive step length in mass-weighted coordinate units.
        fixed_atoms: Optional boolean atom mask; masked coordinates stay fixed.
        initial_direction_q: Optional direction used at the current point,
            useful for choosing either branch from an exact saddle. This is
            not inferred from a Hessian.
        minimum_step_size: Smallest trial size; rejection returns it unchanged.

    Returns:
        Accepted coordinates/energy, or the original coordinates/energy if
        rejected. The returned next step is the attempted size when accepted
        and half that size when rejected.

    Raises:
        ValueError: If inputs, evaluator outputs, or trial step are invalid.
    """
    x0 = np.asarray(coordinates, dtype=float)
    mass = np.asarray(masses, dtype=float)
    if x0.ndim != 2 or x0.shape[1] != 3 or mass.shape != (x0.shape[0],):
        raise ValueError("coordinates must be (n_atoms, 3) and masses (n_atoms,)")
    _finite_array(x0, "coordinates")
    _finite_array(mass, "masses")
    if np.any(mass <= 0.0):
        raise ValueError("masses must be positive")
    if not np.isfinite(step_size) or step_size <= 0.0:
        raise ValueError("step_size must be finite and positive")
    if not np.isfinite(minimum_step_size) or minimum_step_size <= 0.0:
        raise ValueError("minimum_step_size must be finite and positive")

    fixed = (
        np.zeros(x0.shape[0], dtype=bool)
        if fixed_atoms is None
        else np.asarray(fixed_atoms, dtype=bool)
    )
    if fixed.shape != (x0.shape[0],):
        raise ValueError("fixed_atoms must have shape (n_atoms,)")
    energy0, grad0 = energy_gradient(x0.copy())
    grad0 = np.asarray(grad0, dtype=float)
    if not np.isfinite(energy0) or grad0.shape != x0.shape:
        raise ValueError("evaluator must return finite energy and matching gradient")
    _finite_array(grad0, "gradient")
    if initial_direction_q is not None:
        initial = np.asarray(initial_direction_q, dtype=float)
        if initial.shape != x0.shape:
            raise ValueError("initial_direction_q must match coordinates")
        _finite_array(initial, "initial_direction_q")
        initial = initial.copy()
        initial[fixed] = 0.0
        norm = float(np.linalg.norm(initial))
        if norm == 0.0:
            raise ValueError("initial_direction_q has no movable component")
        direction0 = initial / norm
    else:
        direction0 = _direction(grad0, mass, fixed)

    root_mass = np.sqrt(mass)[:, None]
    q0 = root_mass * x0
    q_predict = q0 + step_size * direction0
    x_predict = q_predict / root_mass
    x_predict[fixed] = x0[fixed]
    energy_predict, grad_predict = energy_gradient(x_predict.copy())
    grad_predict = np.asarray(grad_predict, dtype=float)
    if not np.isfinite(energy_predict) or grad_predict.shape != x0.shape:
        raise ValueError("evaluator returned invalid predictor energy or gradient")
    _finite_array(grad_predict, "predictor gradient")
    direction_predict = _direction(grad_predict, mass, fixed)

    q_new = q0 + 0.5 * step_size * (direction0 + direction_predict)
    x_new = q_new / root_mass
    x_new[fixed] = x0[fixed]
    energy_new, grad_new = energy_gradient(x_new.copy())
    grad_new = np.asarray(grad_new, dtype=float)
    if not np.isfinite(energy_new) or grad_new.shape != x0.shape:
        raise ValueError("evaluator returned invalid corrected energy or gradient")
    _finite_array(grad_new, "corrected gradient")
    if energy_new < energy0:
        return DescentStep(x_new, float(energy_new), float(step_size), True)
    return DescentStep(
        x0.copy(), float(energy0), max(step_size / 2.0, minimum_step_size), False
    )
