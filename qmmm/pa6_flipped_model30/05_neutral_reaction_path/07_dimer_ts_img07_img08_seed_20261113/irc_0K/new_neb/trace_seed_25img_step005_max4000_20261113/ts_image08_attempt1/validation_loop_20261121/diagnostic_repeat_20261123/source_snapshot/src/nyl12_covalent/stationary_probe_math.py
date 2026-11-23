"""Numerical summaries for independent stationary-candidate curvature probes."""

from __future__ import annotations

import numpy as np


def curvature_record(
    plus_gradient: np.ndarray,
    minus_gradient: np.ndarray,
    direction: np.ndarray,
    step: float,
    plus_energy: float,
    minus_energy: float,
    center_energy: float,
    expected: float,
) -> dict[str, float]:
    """Compare central directional gradient and energy curvatures.

    Args:
        plus_gradient: Cartesian gradient at positive displacement.
        minus_gradient: Cartesian gradient at negative displacement.
        direction: Cartesian displacement per mass-weighted mode coordinate.
        step: Positive mass-weighted displacement amplitude.
        plus_energy: Positive-displacement energy in Hartree.
        minus_energy: Negative-displacement energy in Hartree.
        center_energy: Undisplaced energy in Hartree.
        expected: Saved mass-weighted Hessian eigenvalue.

    Returns:
        Independent curvature estimates and branch energy differences.

    Raises:
        ValueError: If arrays, numbers or displacement are invalid.
    """
    arrays = [np.asarray(x) for x in (plus_gradient, minus_gradient, direction)]
    numbers = np.array([step, plus_energy, minus_energy, center_energy, expected])
    if (
        step <= 0
        or not np.isfinite(numbers).all()
        or len({x.shape for x in arrays}) != 1
        or not all(np.isfinite(x).all() for x in arrays)
    ):
        raise ValueError("Invalid finite-difference probe")
    gradient_curvature = float(np.sum((arrays[0] - arrays[1]) * arrays[2]) / (2 * step))
    return {
        "step_mass_weighted": step,
        "gradient_curvature": gradient_curvature,
        "energy_curvature": (plus_energy + minus_energy - 2 * center_energy) / step**2,
        "saved_eigenvalue": expected,
        "gradient_curvature_difference": gradient_curvature - expected,
        "plus_energy_difference_hartree": plus_energy - center_energy,
        "minus_energy_difference_hartree": minus_energy - center_energy,
    }
