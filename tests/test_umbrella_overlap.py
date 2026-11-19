"""Tests for umbrella-coordinate parsing and overlap diagnostics."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.umbrella_overlap import (
    campaign_samples,
    overlap_summary,
    read_dumpave,
    write_overlap_data,
)


def _campaign(tmp_path: Path, complete: bool = True) -> Path:
    windows = [
        {"label": "window_00", "target_A": -0.5},
        {"label": "window_01", "target_A": 0.0},
        {"label": "window_02", "target_A": 0.5},
    ]
    (tmp_path / "inputs.json").write_text(json.dumps({"windows": windows}))
    for index, window in enumerate(windows if complete else windows[:-1]):
        stage = tmp_path / window["label"] / "03_production_2ps"
        stage.mkdir(parents=True)
        values = [index * 0.4 - 0.2, index * 0.4, index * 0.4 + 0.2]
        (stage / "umbrella_coordinate.dat").write_text(
            "".join(
                f"{step * 10:8d} {value:9.3f}\n" for step, value in enumerate(values)
            )
        )
        (stage / "process.json").write_text(json.dumps({"exit_status": 0}))
    return tmp_path


def test_dumpave_and_neighbor_overlap(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    source = campaign / "window_00/03_production_2ps/umbrella_coordinate.dat"
    assert read_dumpave(source).tolist() == [-0.2, 0.0, 0.2]
    windows = campaign_samples(campaign)
    summaries, edges = overlap_summary(windows, bins=6)
    assert len(summaries) == 3
    assert summaries[0]["neighbor_overlap"] > 0.0
    assert summaries[-1]["neighbor_overlap"] is None
    output = tmp_path / "analysis"
    output.mkdir()
    json_path, csv_path = write_overlap_data(summaries, edges, output)
    assert json.loads(json_path.read_text())["windows"][0]["sample_count"] == 3
    assert csv_path.read_text().startswith("window,")


def test_incomplete_campaign_fails_unless_explicit(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path, complete=False)
    with pytest.raises(ValueError, match="Incomplete"):
        campaign_samples(campaign)
    assert len(campaign_samples(campaign, allow_incomplete=True)) == 2


def test_discard_fraction_validation(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    with pytest.raises(ValueError, match="discard_fraction"):
        campaign_samples(campaign, discard_fraction=1.0)


def test_campaign_samples_concatenates_completed_production_blocks(
    tmp_path: Path,
) -> None:
    campaign = _campaign(tmp_path)
    for window in ("window_00", "window_01", "window_02"):
        block = campaign / window / "04_production_block_02_5ps"
        block.mkdir()
        (block / "umbrella_coordinate.dat").write_text("0 1.0\n10 1.1\n")
        (block / "process.json").write_text(json.dumps({"exit_status": 0}))
    rows = campaign_samples(campaign, stage="*production*")
    assert all(row["samples_A"].size == 5 for row in rows)
    assert all(";" in row["source"] for row in rows)
