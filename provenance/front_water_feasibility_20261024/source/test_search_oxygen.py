"""Focused tests for physically bounded front-side oxygen feasibility."""

from pathlib import Path
import json
import shutil

import diagnose_oxygen
import numpy as np
import pytest
import search_oxygen

from search_oxygen import front_frame, oxygen_grid, screen_oxygen


def test_grid_distances_front_and_rigid_transform() -> None:
    n, og, preferred = (
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 2.0]),
        np.array([1.0, 0.0, 1.0]),
    )
    xyz, front = oxygen_grid(n, og, preferred)
    assert len(xyz) == 225000
    for center in (n, og):
        distance = np.linalg.norm(xyz - center, axis=1)
        assert distance.min() >= 2.6 - 1e-10
        assert distance.max() <= 3.2 + 1e-10
    assert np.all(front > 0)
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.array([2.0, 4.0, -5.0])
    moved, moved_front = oxygen_grid(
        *(p @ rotation.T + translation for p in (n, og, preferred))
    )
    np.testing.assert_allclose(moved, xyz @ rotation.T + translation, atol=1e-12)
    np.testing.assert_allclose(front, moved_front)


def test_undefined_front_fails() -> None:
    with pytest.raises(ValueError, match="front direction"):
        front_frame(np.zeros(3), np.array([0.0, 0.0, 2.0]), np.array([0.0, 0.0, 1.0]))


def test_all_atom_clearance_and_far_atom_bound() -> None:
    environment = (
        np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]]),
        np.array([1.7, 1.8]),
        ["ASP:272:CA", "MET:20:SD"],
        np.full((2, 3), np.nan),
    )
    margins, collider, local = screen_oxygen(
        np.array([[1.44, 0.0, 0.0], [3.0, 0.0, 0.0]]), environment, np.zeros(3)
    )
    np.testing.assert_allclose(margins, [-1.38, 0.18])
    np.testing.assert_array_equal(collider, [0, 0])
    np.testing.assert_array_equal(local, [0])


def test_hydrogen_bond_exception_requires_angle_and_parent_distance() -> None:
    environment = (
        np.array([[0.0, 0.0, 0.0]]),
        np.array([1.2]),
        ["PHE:113:H"],
        np.array([[-1.0, 0.0, 0.0]]),
    )
    margins, _, _ = screen_oxygen(
        np.array([[1.8, 0.0, 0.0], [0.0, 1.8, 0.0], [0.5, 0.0, 0.0]]),
        environment,
        np.zeros(3),
    )
    assert margins[0] == pytest.approx(0.3)
    assert margins[1] == pytest.approx(-0.52)
    assert margins[2] == pytest.approx(-1.82)


def test_intended_heavy_atom_exception_is_specific() -> None:
    points = np.array([[2.6, 0.0, 0.0]])
    environment = (
        np.zeros((2, 3)),
        np.array([1.55, 1.55]),
        ["T1P:233:N", "ASN:186:ND2"],
        np.full((2, 3), np.nan),
    )
    margins, collider, _ = screen_oxygen(points, environment, np.zeros(3))
    assert margins[0] == pytest.approx(-0.07)
    assert collider[0] == 1


def test_real_source_mapping_and_report(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(search_oxygen, "OUT", tmp_path)
    search_oxygen.main()
    report = json.loads((tmp_path / "oxygen_feasibility.json").read_text())
    assert report["grid_count"] == 225000
    assert report["environment_count"] > 5000
    assert report["unchanged_nonwater_heavy_coordinates"]
    assert len(report["source_sha256"]) == 8


def test_cutoff_cannot_silently_omit_atoms() -> None:
    env = (np.zeros((1, 3)), np.array([1.7]), ["ASP:272:CA"], np.full((1, 3), np.nan))
    with pytest.raises(ValueError, match="outside declared screening bound"):
        screen_oxygen(np.array([[3.3, 0.0, 0.0]]), env, np.zeros(3))


def test_diagnostics_preserve_failed_acceptance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    shutil.copy2(search_oxygen.OUT / "oxygen_grid.npz", tmp_path / "oxygen_grid.npz")
    monkeypatch.setattr(diagnose_oxygen, "OUT", tmp_path)
    diagnose_oxygen.main()
    report = json.loads((tmp_path / "oxygen_diagnostics.json").read_text())
    assert not report["heavy_only_acceptance"]
    assert not report["widened_3p3_diagnostic"]["accepted"]
    assert report["widened_3p3_diagnostic"]["point_count"] == 302760
    assert len(report["user_oxygen_contacts"]) == 20
