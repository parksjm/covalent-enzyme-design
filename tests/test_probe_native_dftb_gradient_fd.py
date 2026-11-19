"""Tests for active-mode finite-difference displacement preparation."""

from pathlib import Path

import numpy as np
import pytest

from scripts.probe_native_dftb_gradient_fd import displace_active, read_final_energy


def test_displace_active_normalizes_mode_and_preserves_frozen_atoms() -> None:
    coordinates = np.zeros((4, 3))
    active = np.array([1, 3])
    mode = np.array([[3.0, 0.0, 0.0], [0.0, 4.0, 0.0]])

    displaced = displace_active(coordinates, active, mode, 0.1)

    assert np.linalg.norm(displaced[active] - coordinates[active]) == pytest.approx(0.1)
    assert np.array_equal(displaced[[0, 2]], coordinates[[0, 2]])


def test_displace_active_rejects_invalid_mode_shape() -> None:
    with pytest.raises(ValueError, match="match the active atom count"):
        displace_active(np.zeros((2, 3)), np.array([0, 1]), np.zeros((1, 3)), 0.1)


def test_displace_active_rejects_zero_mode() -> None:
    with pytest.raises(ValueError, match="nonzero"):
        displace_active(np.zeros((1, 3)), np.array([0]), np.zeros((1, 3)), 0.1)


def test_read_final_energy_accepts_sander_single_point_label(tmp_path: Path) -> None:
    mdout = tmp_path / "mdout"
    mdout.write_text("Energy calculation finished, energy: -7.18E+02\n")

    assert read_final_energy(mdout) == pytest.approx(-718.0)
