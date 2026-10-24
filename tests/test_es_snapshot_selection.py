"""Validation and deterministic physical-time ES selection tests."""

import csv
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent.es_snapshot_selection import METRIC_COLUMNS, select_es_frames


def _write(path: Path, changes: list[dict[str, Any]] | None = None) -> Path:
    rows = []
    for index in range(1, 7):
        row = dict(
            model="model41",
            replica="replica_01",
            trajectory_kind="full",
            frame=index,
            time_ns_estimate=(index - 1) * 0.5,
        )
        row.update(dict(zip(METRIC_COLUMNS, [3.0, 3.5, 3.5, 100.0, -90.0])))
        rows.append(row)
    for index, change in enumerate(changes or []):
        rows[index].update(change)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _select(path: Path, **kwargs: Any) -> list[dict[str, Any]]:
    options: dict[str, Any] = dict(
        model="model41",
        replica="replica_01",
        windows_ns=[(0.0, 1.0), (1.0, 3.0)],
        full_frame_stride=10,
        full_frame_dt_ns=0.05,
    )
    options.update(kwargs)
    return select_es_frames(path, **options)


def test_physical_index_offset_and_midpoint(tmp_path: Path) -> None:
    result = _select(_write(tmp_path / "geometry.csv"))
    assert [row["ordinal_frame"] for row in result] == [2, 5]
    assert [row["full_frame"] for row in result] == [11, 41]
    assert result[0]["physical_time_ns"] == pytest.approx(0.55)
    assert result[0]["logged_time_ns"] == 0.5
    assert result[0]["time_offset_ns"] == pytest.approx(0.05)
    assert result[0]["total_count"] == result[0]["eligible_count"] == 2


def test_geometry_filters_without_ranking_by_geometry(tmp_path: Path) -> None:
    path = _write(tmp_path / "geometry.csv", [{METRIC_COLUMNS[0]: 1.0}])
    assert _select(path)[0]["frame"] == 2
    path = _write(path, [{}, {METRIC_COLUMNS[0]: 3.51}])
    result = _select(path)
    assert result[0]["frame"] == 1
    assert result[0]["eligible_count"] == 1
    assert result[0]["total_count"] == 2


def test_half_open_windows_and_tie_break(tmp_path: Path) -> None:
    path = _write(tmp_path / "geometry.csv")
    result = _select(path, windows_ns=[(0.05, 0.55), (0.55, 1.55)])
    assert [row["frame"] for row in result] == [1, 3]
    tied = _select(
        path, windows_ns=[(0.5, 2.0)], full_frame_stride=1, full_frame_dt_ns=0.5
    )
    assert tied[0]["frame"] == 2


@pytest.mark.parametrize("value", [float("nan"), float("inf"), "bad", ""])
def test_invalid_numeric_rows(tmp_path: Path, value: Any) -> None:
    path = _write(tmp_path / "geometry.csv", [{METRIC_COLUMNS[0]: value}])
    with pytest.raises(ValueError):
        _select(path)


@pytest.mark.parametrize("frame", [0, -1, 1.5, 2, float("nan")])
def test_invalid_frame_indices(tmp_path: Path, frame: float) -> None:
    path = _write(tmp_path / "geometry.csv", [{"frame": frame}])
    with pytest.raises(ValueError):
        _select(path)


def test_duplicate_and_missing_frames(tmp_path: Path) -> None:
    for frame in (1, 3):
        path = _write(tmp_path / "geometry.csv", [{}, {"frame": frame}])
        with pytest.raises(ValueError, match="sequential"):
            _select(path)


@pytest.mark.parametrize(
    "options",
    [
        {"windows_ns": []},
        {"windows_ns": [(1.0, 0.0)]},
        {"windows_ns": [(-1.0, 1.0)]},
        {"windows_ns": [(0.0, float("nan"))]},
        {"windows_ns": [(0.0, 2.0), (1.0, 3.0)]},
        {"windows_ns": [(10.0, 11.0)]},
        {"full_frame_stride": 0},
        {"full_frame_stride": True},
        {"full_frame_stride": 1.5},
        {"full_frame_dt_ns": 0.0},
        {"model": "other"},
        {"replica": "other"},
    ],
)
def test_invalid_arguments(tmp_path: Path, options: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        _select(_write(tmp_path / "geometry.csv"), **options)


@pytest.mark.parametrize(
    "metric,value",
    [(0, 3.6), (1, 4.1), (2, 4.1), (3, 84.0), (3, 126.0), (4, 59.0), (4, -126.0)],
)
def test_no_eligible_geometry(tmp_path: Path, metric: int, value: float) -> None:
    path = _write(tmp_path / "geometry.csv", [{METRIC_COLUMNS[metric]: value}])
    with pytest.raises(ValueError, match="No eligible"):
        _select(path, windows_ns=[(0.0, 0.5)])


def test_required_columns_and_negative_distance(tmp_path: Path) -> None:
    path = tmp_path / "geometry.csv"
    path.write_text("model,frame\nmodel41,1\n")
    with pytest.raises(ValueError, match="columns"):
        _select(path)
    _write(path, [{METRIC_COLUMNS[0]: -1}])
    with pytest.raises(ValueError, match="positive"):
        _select(path)


def test_ignores_nonmatching_trajectory_kind(tmp_path: Path) -> None:
    path = _write(tmp_path / "geometry.csv", [{"trajectory_kind": "compact"}] * 6)
    with pytest.raises(ValueError, match="No matching"):
        _select(path)
