"""Tests for ordered umbrella-window VMD frame selection."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.io import netcdf_file

from nyl12_covalent.umbrella_vmd import (
    campaign_frame_plan,
    completed_production_trajectories,
    cpptraj_trajectory_input,
    evenly_spaced_frames,
    netcdf_frame_count,
)


def _trajectory(path: Path, frames: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netcdf_file(path, "w") as handle:
        handle.createDimension("frame", frames)
        handle.createDimension("atom", 2)
        handle.createDimension("spatial", 3)
        coordinates = handle.createVariable(
            "coordinates", "f", ("frame", "atom", "spatial")
        )
        coordinates[:] = np.zeros((frames, 2, 3), dtype=np.float32)


def _completed_stage(window: Path, name: str, frames: int) -> Path:
    stage = window / name
    trajectory = stage / "trajectory.nc"
    _trajectory(trajectory, frames)
    (stage / "process.json").write_text(json.dumps({"exit_status": 0}))
    return trajectory


def test_evenly_spaced_frames_respects_limit_and_endpoints() -> None:
    assert evenly_spaced_frames(5, 3) == [1, 3, 5]
    assert evenly_spaced_frames(3, 10) == [1, 2, 3]
    with pytest.raises(ValueError, match="positive"):
        evenly_spaced_frames(5, 0)


def test_completed_production_trajectories_excludes_failed_stage(
    tmp_path: Path,
) -> None:
    window = tmp_path / "window_00"
    first = _completed_stage(window, "02_production_block_01_5ps", 2)
    failed = window / "02_production_block_02_5ps"
    _trajectory(failed / "trajectory.nc", 2)
    (failed / "process.json").write_text(json.dumps({"exit_status": 1}))
    assert completed_production_trajectories(window) == [first]
    assert netcdf_frame_count(first) == 2


def test_campaign_frame_plan_orders_windows_then_downsampled_frames(
    tmp_path: Path,
) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    windows = [
        {"label": "window_00", "target_A": -0.5},
        {"label": "window_01", "target_A": 0.5},
    ]
    (campaign / "inputs.json").write_text(json.dumps({"windows": windows}))
    for row in windows:
        window = campaign / row["label"]
        _completed_stage(window, "02_production_block_01_5ps", 2)
        _completed_stage(window, "02_production_block_02_5ps", 3)
    plan = campaign_frame_plan(campaign, maximum_frames_per_window=3)
    assert [row["window"] for row in plan] == ["window_00"] * 3 + ["window_01"] * 3
    assert [row["window_global_frame"] for row in plan[:3]] == [1, 3, 5]
    assert [row["source_frame"] for row in plan[:3]] == [1, 1, 3]
    text = cpptraj_trajectory_input(
        tmp_path / "system.parm7", plan, tmp_path / "movie.dcd", ":WAT"
    )
    assert text.count("trajin ") == 6
    assert "strip :WAT" in text
    assert text.index("window_00") < text.index("window_01")


def test_campaign_frame_plan_requires_complete_windows(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "inputs.json").write_text(
        json.dumps(
            {
                "windows": [
                    {"label": "window_00", "target_A": 0.0},
                    {"label": "window_01", "target_A": 0.1},
                ]
            }
        )
    )
    _completed_stage(campaign / "window_00", "03_production_2ps", 2)
    with pytest.raises(ValueError, match="window_01"):
        campaign_frame_plan(campaign, 10)
    plan = campaign_frame_plan(campaign, 10, allow_incomplete=True)
    assert {row["window"] for row in plan} == {"window_00"}
