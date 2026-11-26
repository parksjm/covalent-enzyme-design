"""Strict conditional continuation gates for the bounded hydration followon."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def read(path: Path) -> dict[str, Any]:
    """Read a required JSON object."""
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected object: {path}")
    return value


def digest(path: Path) -> str:
    """Return the SHA256 digest of a required file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_hashes(hashes: dict[str, str]) -> None:
    """Reject missing files, empty manifests and changed file identities."""
    if not hashes:
        raise ValueError("Empty source manifest")
    for filename, expected in hashes.items():
        if digest(Path(filename)) != expected:
            raise ValueError(f"Identity changed: {filename}")


def accepted(base: Path, index: int) -> None:
    """Require actual flat native gates, numeric values and restart lineage."""
    config_name = (
        "extension_20ps_config.json" if index <= 4 else "followon_runtime_config.json"
    )
    config = read(base / config_name)
    block = base / f"extend_{index:02d}"
    for name, expected in config["segment_input_sha256"][block.name].items():
        if digest(block / "inputs" / name) != expected:
            raise ValueError(f"Native input identity changed: {block.name}/{name}")
    status = read(block / "status.json")
    process = read(block / "process.json")
    record = read(block / "scientific_acceptance.json")
    numeric, boundary = record["numeric"], record["boundary"]
    temperature = float(numeric["final_temperature_K"])
    frozen = float(boundary["frozen_max_displacement_A"])
    distance = float(record["C_OG_A"])
    previous = base / ("equil" if index == 1 else f"extend_{index - 1:02d}")
    start_hash = digest(previous / "restart.ncrst")
    flags = [
        status["status"] == "accepted",
        status["exit_code"] == 0,
        process["status"] == "process_completed",
        process["exit_code"] == 0,
        status["restart_sha256"] == digest(block / "restart.ncrst"),
        status["start_restart_sha256"] == start_hash,
        process["start_restart_sha256"] == start_hash,
        digest(block / "inputs/restart.ncrst") == start_hash,
        numeric["completed_steps"] == 10000,
        numeric["electronic_history"]["verified_scc_history"] is True,
        math.isfinite(temperature) and 250 <= temperature <= 350,
        math.isfinite(distance) and 0 < distance < 1.8,
        math.isfinite(frozen) and 0 <= frozen <= 1e-6,
        boundary["box_unchanged"] is True,
    ]
    flags.extend(
        record[key] is True
        for key in (
            "scc_verified",
            "native_process_passed",
            "finite_temperature_passed",
            "acyl_retained",
        )
    )
    if not all(flags):
        raise ValueError(f"Native acceptance/provenance failed: {block.name}")
    trajectory = block / "trajectory.nc"
    if not trajectory.is_file() or trajectory.stat().st_size == 0:
        raise ValueError(f"Missing trajectory: {block.name}")


def report_candidates(base: Path, count: int) -> list[dict[str, Any]]:
    """Verify full accepted-trajectory report and return all near-attack waters."""
    analysis = base / "analysis"
    manifest = read(analysis / "source_manifest.json")["sha256"]
    verify_hashes(manifest)
    report = read(analysis / "report.json")
    expected = {f"extend_{index:02d}" for index in range(1, count + 1)}
    stages = report["frame_count_by_stage"]
    observed = {name for name in stages if name.startswith("extend_")}
    if observed != expected or any(stages[name] != 100 for name in expected):
        raise ValueError("Analysis does not cover every accepted 5 ps block")
    sources = {entry["stage"]: entry["path"] for entry in report["trajectory_sources"]}
    for name in expected:
        directory = base / name
        required = [
            directory / "trajectory.nc",
            directory / "scientific_acceptance.json",
        ]
        if sources[name] != str(directory / "trajectory.nc"):
            raise ValueError("Wrong trajectory source")
        if any(str(path) not in manifest for path in required):
            raise ValueError("Report lacks required native source hashes")
    for filename in ("report.json", "per_water_summary.json", "per_frame_metrics.csv"):
        if str(analysis / filename) not in manifest:
            raise ValueError("Report artifact missing from manifest")
    criteria = report["criteria"]
    if criteria["near_attack_O_C_A_max"] != 3.5 or criteria["attack_angle_deg"] != [
        95,
        125,
    ]:
        raise ValueError("Changed candidate geometry criteria")
    with (analysis / "nearest_water_per_frame.csv").open() as stream:
        final = list(csv.DictReader(stream))[-1]
    if not math.isclose(
        float(final["time_ps"]) - report["endpoint_time_ps"],
        6 + count * 5,
        abs_tol=1e-6,
    ):
        raise ValueError(
            "Analysis final time does not cover required hydration duration"
        )
    waters = read(analysis / "per_water_summary.json")["waters"]
    if len(waters) != report["water_count"]:
        raise ValueError("Incomplete per-water summary")
    return [water for water in waters if water["near_attack_frame_count"] > 0]
