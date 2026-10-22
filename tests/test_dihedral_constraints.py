"""Independent torsion derivatives and rigid-motion invariance checks."""

import numpy as np
import pytest

from nyl12_covalent.dihedral_constraints import (
    dihedral_constraint_derivatives,
    rigid_cartesian_directions,
)


def _coordinates(angle: float) -> np.ndarray:
    """Return a torsion whose last atom gives phi=atan2(z,y)."""
    return np.asarray(
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, np.cos(angle), np.sin(angle)],
        ]
    )


@pytest.mark.parametrize("angle", [np.pi / 4, np.pi - 1e-6, -np.pi + 1e-6])
def test_matches_analytic_atan2_derivatives(angle: float) -> None:
    xyz = _coordinates(angle)
    jacobian, hessians = dihedral_constraint_derivatives(
        xyz, [(0, 1, 2, 3)], step_bohr=1e-4
    )
    y, z = xyz[3, 1:]
    np.testing.assert_allclose(jacobian[0, 10:12], [-z, y], atol=2e-8)
    expected = [[2 * y * z, z * z - y * y], [z * z - y * y, -2 * y * z]]
    np.testing.assert_allclose(hessians[0, 10:12, 10:12], expected, atol=2e-7)


def test_derivatives_scale_with_length_and_ignore_rigid_motions() -> None:
    xyz = _coordinates(0.7)
    j1, h1 = dihedral_constraint_derivatives(xyz, [(0, 1, 2, 3)], step_bohr=2e-4)
    j2, h2 = dihedral_constraint_derivatives(2 * xyz, [(0, 1, 2, 3)], step_bohr=2e-4)
    np.testing.assert_allclose(j1, 2 * j2, atol=1e-7)
    np.testing.assert_allclose(h1, 4 * h2, atol=4e-7)
    rigid = rigid_cartesian_directions(xyz)
    assert rigid.shape == (12, 6)
    np.testing.assert_allclose(j1 @ rigid, 0, atol=1e-7)
    np.testing.assert_allclose(h1[0], h1[0].T, atol=0)


@pytest.mark.parametrize("step", [0.0, -1.0, np.nan])
def test_rejects_invalid_steps(step: float) -> None:
    with pytest.raises(ValueError, match="step"):
        dihedral_constraint_derivatives(
            _coordinates(0.5), [(0, 1, 2, 3)], step_bohr=step
        )


@pytest.mark.parametrize(
    "torsions", [[], [(0, 1, 2, 3)] * 2, [(0, 0, 2, 3)], [(0, 1, 2, 4)]]
)
def test_rejects_invalid_torsions(torsions: list[tuple[int, int, int, int]]) -> None:
    with pytest.raises(ValueError):
        dihedral_constraint_derivatives(_coordinates(0.5), torsions, step_bohr=1e-4)


def test_rejects_collinear_torsion_and_linear_rigid_basis() -> None:
    xyz = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    )
    with pytest.raises(ValueError, match="Collinear"):
        dihedral_constraint_derivatives(xyz, [(0, 1, 2, 3)], step_bohr=1e-4)
    with pytest.raises(ValueError, match="nonlinear"):
        rigid_cartesian_directions(xyz)


def test_rejects_nonfinite_coordinates() -> None:
    xyz = _coordinates(0.5)
    xyz[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        dihedral_constraint_derivatives(xyz, [(0, 1, 2, 3)], step_bohr=1e-4)
    with pytest.raises(ValueError, match="finite"):
        rigid_cartesian_directions(xyz)
