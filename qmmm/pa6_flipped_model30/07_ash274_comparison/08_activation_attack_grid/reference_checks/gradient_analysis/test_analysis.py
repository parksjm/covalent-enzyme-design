"""Focused saved-gradient parser and projection checks."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "saved_analysis", Path(__file__).with_name("analyze.py")
)
a = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(a)


def test_dual_directions() -> None:
    """Coupled CVs obtain unit changes without cross-coordinate motion."""
    j = np.array([[[1.0, 2.0, 0.0]], [[2.0, 0.0, 1.0]]])
    d = a.dual_directions(j)
    np.testing.assert_allclose(
        j.reshape(2, -1) @ d.reshape(2, -1).T, np.eye(2), atol=1e-14
    )
    with pytest.raises(ValueError):
        a.dual_directions(np.array([j[0], j[0]]))


def test_fchk(tmp_path: Path) -> None:
    """Declared arrays parse Fortran exponents and reject missing labels."""
    p = tmp_path / "test.fchk"
    p.write_text("Cartesian Gradient R N= 3\n1D-2 2D-2 3D-2\n")
    np.testing.assert_allclose(
        a.fchk_array(p, "Cartesian Gradient"), [0.01, 0.02, 0.03]
    )
    with pytest.raises(ValueError):
        a.fchk_array(p, "Absent")


def test_engrad_units(tmp_path: Path) -> None:
    """xTB coordinates are bohr; gradients remain Eh/bohr until conversion."""
    p = tmp_path / "test.engrad"
    p.write_text("# comment\n1\n-1.0\n0.1\n0.2\n0.3\n1 1 0 0\n")
    z, x, g = a.engrad(p)
    assert z.tolist() == [1]
    assert x[0, 0] == a.BOHR
    np.testing.assert_allclose(g, [[0.1, 0.2, 0.3]])
    assert a.HARTREE / a.BOHR == pytest.approx(1185.8210466, rel=1e-8)


def test_distance() -> None:
    """Distance Jacobian has correct signed derivative."""
    r, j = a.distance(np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]), 0, 1)
    assert r == 2
    np.testing.assert_allclose(j, [[-1, 0, 0], [1, 0, 0]])
