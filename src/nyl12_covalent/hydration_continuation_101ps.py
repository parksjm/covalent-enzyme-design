"""Acceptance and analysis gates for the 101 ps hydration continuation."""

from __future__ import annotations

import hashlib
import json
import math

import numpy as np
from pathlib import Path


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_accepted_segment(base: Path, index: int) -> None:
    """Validate a contiguous native segment, input hashes, chemistry and restart."""
    config_name = (
        "extension_20ps_config.json"
        if index <= 4
        else "followon_runtime_config.json"
        if index <= 9
        else "continuation_101ps_config.json"
    )
    config = json.loads((base / config_name).read_text())
    name = f"extend_{index:02d}"
    block = base / name
    status = json.loads((block / "status.json").read_text())
    process = json.loads((block / "process.json").read_text())
    acceptance = json.loads((block / "scientific_acceptance.json").read_text())
    numeric, boundary = acceptance["numeric"], acceptance["boundary"]
    expected_start = (
        base / ("equil" if index == 1 else f"extend_{index - 1:02d}") / "restart.ncrst"
    )
    start_hash = sha256(expected_start)
    conditions = [
        status.get("status") == "accepted",
        status.get("exit_code") == 0,
        process.get("status") == "process_completed",
        process.get("exit_code") == 0,
        numeric.get("completed_steps") == 10000,
        numeric.get("electronic_history", {}).get("verified_scc_history") is True,
        250 <= float(numeric["final_temperature_K"]) <= 350,
        0 < float(acceptance["C_OG_A"]) < 1.8,
        acceptance.get("scc_verified") is True,
        acceptance.get("acyl_retained") is True,
        acceptance.get("finite_temperature_passed") is True,
        acceptance.get("native_process_passed") is True,
        boundary.get("box_unchanged") is True,
        float(boundary.get("frozen_max_displacement_A", 1)) <= 1e-6,
        status.get("restart_sha256") == sha256(block / "restart.ncrst"),
        status.get("start_restart_sha256") == start_hash,
        process.get("start_restart_sha256") == start_hash,
        sha256(block / "inputs/restart.ncrst") == start_hash,
    ]
    for filename, expected in config["segment_input_sha256"][name].items():
        conditions.append(sha256(block / "inputs" / filename) == expected)
    if not all(conditions):
        raise ValueError(f"Native acceptance or restart lineage failed: {name}")
    trajectory = block / "trajectory.nc"
    if not trajectory.is_file() or trajectory.stat().st_size == 0:
        raise ValueError(f"Missing accepted trajectory: {name}")


def validate_prefix(base: Path, first: int, last: int) -> int:
    """Require all historical blocks and only a contiguous attempted prefix."""
    for index in range(1, first):
        validate_accepted_segment(base, index)
    accepted_count = 0
    gap = False
    for index in range(first, last + 1):
        marker = base / f"extend_{index:02d}/status.json"
        if not marker.exists():
            gap = True
            continue
        status = json.loads(marker.read_text())
        if status.get("status") == "accepted":
            if gap:
                raise ValueError("Accepted continuation blocks are not contiguous")
            validate_accepted_segment(base, index)
            accepted_count += 1
        else:
            gap = True
    return accepted_count


def candidate_found(base: Path, count: int) -> bool:
    """Require fresh complete analysis and screen the reviewed near-attack gate."""
    manifest = json.loads((base / "analysis/source_manifest.json").read_text())[
        "sha256"
    ]
    for filename, expected in manifest.items():
        if sha256(Path(filename)) != expected:
            raise ValueError(f"Analysis source changed: {filename}")
    report = json.loads((base / "analysis/report.json").read_text())
    stages = report["frame_count_by_stage"]
    expected = {f"extend_{i:02d}" for i in range(1, count + 1)}
    observed = {name for name in stages if name.startswith("extend_")}
    if observed != expected or any(stages[name] != 100 for name in expected):
        raise ValueError("Analysis does not cover the current accepted segment count")
    criteria = report["criteria"]
    if criteria["near_attack_O_C_A_max"] != 3.5 or criteria["attack_angle_deg"] != [
        95,
        125,
    ]:
        raise ValueError("Near-attack screen changed")
    waters = json.loads((base / "analysis/per_water_summary.json").read_text())[
        "waters"
    ]
    if len(waters) != report["water_count"]:
        raise ValueError("Per-water analysis is incomplete")
    return any(item["near_attack_frame_count"] > 0 for item in waters)


def native_time_matches(observed_ps: float, expected_ps: float) -> bool:
    """Compare serialized NetCDF float32 time with its scheduled endpoint."""
    scale = max(abs(observed_ps), abs(expected_ps), 1.0)
    ulp = abs(float(np.spacing(np.float32(scale))))
    tolerance = max(
        1e-6,
        ulp / 2 + 8 * np.finfo(np.float64).eps * scale,
    )
    return math.isclose(observed_ps, expected_ps, rel_tol=0.0, abs_tol=tolerance)
