"""Read accepted native hydration blocks without modifying their receipts."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def accepted_extensions(base: Path) -> list[Path]:
    """Validate the native extension schema and return its contiguous prefix.

    Args:
        base: Hydration directory containing extend_01 through extend_04.

    Returns:
        Ordered scientifically accepted segment directories.

    Raises:
        ValueError: If a claimed accepted block fails the original native gates.
    """
    accepted = []
    pending_seen = False
    for index in range(1, 5):
        segment = base / f"extend_{index:02d}"
        marker = segment / "status.json"
        if not marker.exists():
            pending_seen = True
            continue
        status = json.loads(marker.read_text())
        if status.get("status") != "accepted":
            pending_seen = True
            continue
        if pending_seen:
            raise ValueError("Accepted hydration blocks are not a contiguous prefix")
        record = json.loads((segment / "scientific_acceptance.json").read_text())
        process = json.loads((segment / "process.json").read_text())
        numeric, boundary = record["numeric"], record["boundary"]
        temperature = float(numeric["final_temperature_K"])
        distance = float(record["C_OG_A"])
        frozen = float(boundary["frozen_max_displacement_A"])
        conditions = [
            numeric["completed_steps"] == 10000,
            math.isfinite(temperature) and 250 <= temperature <= 350,
            numeric["electronic_history"]["verified_scc_history"] is True,
            math.isfinite(distance) and 0 < distance < 1.8,
            math.isfinite(frozen) and 0 <= frozen <= 1e-6,
            boundary["box_unchanged"] is True,
            process["status"] == "process_completed",
            process["exit_code"] == 0,
        ]
        conditions.extend(
            record[name] is True
            for name in (
                "scc_verified",
                "acyl_retained",
                "finite_temperature_passed",
                "native_process_passed",
            )
        )
        if not all(conditions):
            raise ValueError(f"Native hydration acceptance failed: {segment.name}")
        restart = segment / "restart.ncrst"
        if hashlib.sha256(restart.read_bytes()).hexdigest() != status["restart_sha256"]:
            raise ValueError(f"Accepted restart identity changed: {segment.name}")
        trajectory = segment / "trajectory.nc"
        if not trajectory.is_file() or trajectory.stat().st_size == 0:
            raise ValueError(f"Accepted trajectory is missing/empty: {segment.name}")
        accepted.append(segment)
    return accepted
