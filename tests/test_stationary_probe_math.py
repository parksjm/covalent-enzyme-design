"""Analytic harmonic checks for the independent curvature summaries."""

import numpy as np
import pytest

from nyl12_covalent.stationary_probe_math import curvature_record


@pytest.mark.parametrize("eigenvalue", [-0.23, 3.4e-5, 0.5])
def test_harmonic_curvatures_with_nonunit_cartesian_mode(eigenvalue: float) -> None:
    """Mass-scaled directions retain the specified Rayleigh curvature."""
    direction = np.array([[0.25, 0.0, 0.0]])
    step = 0.025
    hessian = eigenvalue / 0.25**2
    gradient = hessian * step * direction
    energy = 0.5 * eigenvalue * step**2
    result = curvature_record(
        gradient, -gradient, direction, step, energy, energy, 0.0, eigenvalue
    )
    assert result["gradient_curvature"] == pytest.approx(eigenvalue)
    assert result["energy_curvature"] == pytest.approx(eigenvalue)
    assert result["gradient_curvature_difference"] == pytest.approx(0.0)


@pytest.mark.parametrize("step, bad", [(0.0, False), (-1.0, False), (0.1, True)])
def test_invalid_probe_rejected(step: float, bad: bool) -> None:
    """Invalid displacements and nonfinite gradients fail explicitly."""
    gradient = np.array([np.nan if bad else 0.0])
    with pytest.raises(ValueError):
        curvature_record(gradient, gradient, np.ones(1), step, 0.0, 0.0, 0.0, 0.0)


def test_shape_mismatch_rejected() -> None:
    """Different gradient and direction shapes cannot silently broadcast."""
    with pytest.raises(ValueError):
        curvature_record(np.zeros(2), np.zeros(2), np.ones(1), 0.1, 0.0, 0.0, 0.0, 0.0)
