"""Tests for the modular WT catalytic-Thr torsion screen."""

from pathlib import Path

import pytest

from scripts.screen_modular_wt_thr_torsions import (
    cpptraj_text,
    read_cpptraj_scalar,
    representative_joint_frames,
    wrapped_delta,
)


@pytest.mark.parametrize(
    "value,reference,expected",
    [(170.0, -170.0, -20.0), (-170.0, 170.0, 20.0), (10.0, 20.0, -10.0)],
)
def test_wrapped_delta(value: float, reference: float, expected: float) -> None:
    assert wrapped_delta(value, reference) == pytest.approx(expected)


def test_read_cpptraj_scalar(tmp_path: Path) -> None:
    path = tmp_path / "values.dat"
    path.write_text("#Frame value\n 1 3.0\n 2 4.0\n", encoding="utf-8")
    assert read_cpptraj_scalar(path) == [3.0, 4.0]
    path.write_text("1 3.0\n3 4.0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Nonsequential"):
        read_cpptraj_scalar(path)


def test_cpptraj_text_declares_thr_torsions(tmp_path: Path) -> None:
    text = cpptraj_text(tmp_path / "top", tmp_path / "traj", 2000, tmp_path)
    assert "trajin" in text and "1 2000 1" in text
    assert "dihedral psi :929@N :929@CA :929@C :930@N" in text
    assert "dihedral chi1 :929@N :929@CA :929@CB :929@OG1" in text
    with pytest.raises(ValueError):
        cpptraj_text(tmp_path / "top", tmp_path / "traj", 0, tmp_path)


def test_representative_joint_frames_prefers_tier_a_then_midpoint() -> None:
    rows = [
        {
            "frame": 1,
            "time_ns": 55.0,
            "tier": "B",
            "native_torsions": True,
            "psi_delta": 1.0,
            "chi1_delta": 2.0,
        },
        {
            "frame": 2,
            "time_ns": 56.0,
            "tier": "A",
            "native_torsions": True,
            "psi_delta": 3.0,
            "chi1_delta": 4.0,
        },
        {
            "frame": 3,
            "time_ns": 60.0,
            "tier": "A",
            "native_torsions": True,
            "psi_delta": 5.0,
            "chi1_delta": 6.0,
        },
        {
            "frame": 4,
            "time_ns": 75.0,
            "tier": "A",
            "native_torsions": False,
            "psi_delta": 7.0,
            "chi1_delta": 8.0,
        },
    ]
    result = representative_joint_frames(rows)
    assert [row["frame"] for row in result] == [3]
    assert result[0]["tier"] == "A"
    assert result[0]["joint_eligible_count"] == 2
