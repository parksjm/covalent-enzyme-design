"""Focused opposite-hemisphere diagnostics verification."""

import json
import shutil
from pathlib import Path

import numpy as np
import opposite_oxygen
import pytest


def test_reflection_is_involutory_and_preserves_axial_distances() -> None:
    n = np.array([1.0, 2.0, 3.0])
    og = n + np.array([0.0, 0.0, 2.0])
    axis = (og - n) / np.linalg.norm(og - n)
    points = n + np.array([[1.0, 2.0, 1.0], [3.0, 1.0, -1.0]])
    reflected = opposite_oxygen.reflect_radial(points, n, axis)
    np.testing.assert_allclose(
        opposite_oxygen.reflect_radial(reflected, n, axis), points
    )
    for center in (n, og):
        np.testing.assert_allclose(
            np.linalg.norm(reflected - center, axis=1),
            np.linalg.norm(points - center, axis=1),
        )
    np.testing.assert_allclose(reflected[:, :2] - n[:2], -(points[:, :2] - n[:2]))


def test_real_opposite_screen_remains_diagnostic(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    shutil.copy2(opposite_oxygen.OUT / "oxygen_grid.npz", tmp_path / "oxygen_grid.npz")
    monkeypatch.setattr(opposite_oxygen, "OUT", tmp_path)
    opposite_oxygen.main()
    report = json.loads((tmp_path / "opposite_oxygen_feasibility.json").read_text())
    assert report["grid_count"] == 225000
    assert report["diagnostic_only"]
    assert report["not_mechanism_evidence"]
    assert not report["candidate_seed_emitted"]
    grid = np.load(tmp_path / "opposite_oxygen_grid.npz")
    assert np.all(grid["front"] < 0)
    assert report["full_pass_count"] == np.sum(grid["full_margin"] >= -1e-8)
    assert report["heavy_pass_count"] == np.sum(grid["heavy_margin"] >= -1e-8)
