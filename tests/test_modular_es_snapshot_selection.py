"""Tests for modular PA6 E-S snapshot selection."""

import csv
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent.modular_es_snapshot_selection import (
    classify_geometry,
    select_modular_es_frames,
)


def _row(frame: int, **changes: float) -> dict[str, float]:
    row = {
        "frame": float(frame),
        "time_ns": (frame - 1) * 0.5,
        "d1": 3.0,
        "d2": 3.0,
        "d3": 3.0,
        "attack": 105.0,
        "ligand_rmsd": 2.0,
        "projected_phi": -90.0,
    }
    row.update(changes)
    return row


def _write(path: Path, changes: dict[int, dict[str, Any]] | None = None) -> Path:
    rows = [_row(frame) for frame in range(1, 7)]
    for index, values in (changes or {}).items():
        rows[index].update(values)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_prefers_tier_a_then_time_midpoint(tmp_path: Path) -> None:
    path = _write(tmp_path / "geometry.csv", {2: {"d1": 3.4}})
    result = select_modular_es_frames(path, [(0.5, 2.0)])
    assert result[0]["frame"] == 4
    assert result[0]["tier"] == "A"
    assert result[0]["window_tier_counts"] == {"A": 2, "B": 1}


def test_falls_back_to_tier_b_without_geometry_ranking(tmp_path: Path) -> None:
    changes = {index: {"d1": 3.4} for index in range(6)}
    changes[1]["d1"] = 3.49
    result = select_modular_es_frames(
        _write(tmp_path / "geometry.csv", changes), [(0.0, 1.5)]
    )
    assert result[0]["frame"] == 2
    assert result[0]["tier"] == "B"


@pytest.mark.parametrize(
    "changes",
    [
        {"d1": 3.51},
        {"d2": 4.01},
        {"attack": 84.9},
        {"projected_phi": 125.1},
    ],
)
def test_classification_rejects_outside_tier_b(changes: dict[str, float]) -> None:
    assert classify_geometry(_row(1, **changes)) is None


@pytest.mark.parametrize(
    "windows",
    [[], [(-1.0, 1.0)], [(1.0, 1.0)], [(0.0, 2.0), (1.0, 3.0)]],
)
def test_rejects_invalid_windows(
    tmp_path: Path, windows: list[tuple[float, float]]
) -> None:
    with pytest.raises(ValueError):
        select_modular_es_frames(_write(tmp_path / "geometry.csv"), windows)


def test_rejects_bad_schema_sequence_and_numbers(tmp_path: Path) -> None:
    path = tmp_path / "geometry.csv"
    path.write_text("frame,time_ns\n1,0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing"):
        select_modular_es_frames(path, [(0.0, 1.0)])
    path = _write(path, {1: {"frame": 7}})
    with pytest.raises(ValueError, match="sequential"):
        select_modular_es_frames(path, [(0.0, 1.0)])
    path = _write(path, {0: {"d1": "nan"}})
    with pytest.raises(ValueError, match="Nonfinite"):
        select_modular_es_frames(path, [(0.0, 1.0)])


def test_rejects_window_without_eligible_geometry(tmp_path: Path) -> None:
    path = _write(tmp_path / "geometry.csv", {index: {"d1": 4.0} for index in range(6)})
    with pytest.raises(ValueError, match="No Tier"):
        select_modular_es_frames(path, [(0.0, 2.0)])
