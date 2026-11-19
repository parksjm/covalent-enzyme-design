"""Select ordered, downsampled frames from QM/MM umbrella campaigns."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import netcdf_file


def netcdf_frame_count(path: Path) -> int:
    """Return the number of coordinate frames without loading coordinates.

    Args:
        path: Amber NetCDF trajectory.

    Returns:
        Positive trajectory frame count.
    """
    with netcdf_file(path, "r", mmap=False) as trajectory:
        coordinates = trajectory.variables.get("coordinates")
        if coordinates is None or len(coordinates.shape) != 3:
            raise ValueError(f"Missing coordinate array in {path}")
        frames = int(coordinates.shape[0])
    if frames < 1:
        raise ValueError(f"Trajectory contains no frames: {path}")
    return frames


def evenly_spaced_frames(frame_count: int, maximum_frames: int) -> list[int]:
    """Choose at most ``maximum_frames`` evenly spaced one-based frames.

    Args:
        frame_count: Available frames.
        maximum_frames: Maximum frames retained for one window.

    Returns:
        Strictly increasing one-based frame indices including both endpoints.
    """
    if frame_count < 1 or maximum_frames < 1:
        raise ValueError("Frame counts and limits must be positive")
    retained = min(frame_count, maximum_frames)
    if retained == frame_count:
        return list(range(1, frame_count + 1))
    indices = np.rint(np.linspace(1, frame_count, retained)).astype(int).tolist()
    if len(indices) != len(set(indices)):
        raise ValueError("Even frame selection unexpectedly produced duplicates")
    return indices


def completed_production_trajectories(window: Path) -> list[Path]:
    """Find completed production trajectories for one umbrella window.

    Args:
        window: Umbrella window directory.

    Returns:
        Production trajectories ordered by stage name.
    """
    trajectories: list[Path] = []
    for stage in sorted(window.iterdir() if window.is_dir() else []):
        if not stage.is_dir() or "production" not in stage.name:
            continue
        process = stage / "process.json"
        trajectory = stage / "trajectory.nc"
        if not process.is_file() or not trajectory.is_file():
            continue
        record = json.loads(process.read_text())
        if record.get("exit_status") == 0:
            trajectories.append(trajectory)
    return trajectories


def campaign_frame_plan(
    campaign: Path,
    maximum_frames_per_window: int,
    allow_incomplete: bool = False,
) -> list[dict[str, Any]]:
    """Build an ordered frame plan across umbrella windows and production blocks.

    Args:
        campaign: Umbrella campaign with ``inputs.json`` and window directories.
        maximum_frames_per_window: Downsampling ceiling applied independently.
        allow_incomplete: Skip windows lacking completed production trajectories.

    Returns:
        Ordered records mapping output frames to source window and trajectory.
    """
    if maximum_frames_per_window < 1:
        raise ValueError("maximum_frames_per_window must be positive")
    inputs = json.loads((campaign / "inputs.json").read_text())
    plan: list[dict[str, Any]] = []
    missing: list[str] = []
    for window in inputs["windows"]:
        label = str(window["label"])
        trajectories = completed_production_trajectories(campaign / label)
        if not trajectories:
            if not allow_incomplete:
                missing.append(label)
            continue
        counts = [netcdf_frame_count(path) for path in trajectories]
        selected = evenly_spaced_frames(sum(counts), maximum_frames_per_window)
        offsets = np.cumsum([0, *counts])
        for global_frame in selected:
            segment_index = int(np.searchsorted(offsets[1:], global_frame, side="left"))
            local_frame = global_frame - int(offsets[segment_index])
            plan.append(
                {
                    "output_frame": len(plan),
                    "window": label,
                    "target_A": float(window["target_A"]),
                    "source": trajectories[segment_index],
                    "source_frame": local_frame,
                    "window_global_frame": global_frame,
                    "window_available_frames": int(sum(counts)),
                }
            )
    if missing:
        raise ValueError(
            "Windows lack completed production trajectories: " + ", ".join(missing)
        )
    if not plan:
        raise ValueError("No completed production frames are available")
    return plan


def cpptraj_trajectory_input(
    topology: Path,
    plan: list[dict[str, Any]],
    output: Path,
    strip_mask: str,
) -> str:
    """Return cpptraj input that concatenates the declared frames in order."""
    if not plan:
        raise ValueError("Frame plan is empty")
    lines = [f"parm {topology.resolve()}"]
    lines.extend(
        f"trajin {row['source'].resolve()} {row['source_frame']} {row['source_frame']}"
        for row in plan
    )
    lines.extend(
        [
            f"strip {strip_mask}",
            f"trajout {output.resolve()} dcd",
            "run",
            "quit",
        ]
    )
    return "\n".join(lines) + "\n"
