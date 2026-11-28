"""Wait for matched WT completion, then run one explicitly frozen workflow."""

from __future__ import annotations

import json
import logging
import math
import os
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from nyl12_covalent.matched_wt_gate import (
    ReceiptReader,
    inspect_cases,
    verify_sources,
    verify_topologies,
)

LOGGER = logging.getLogger(__name__)


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write(path: Path, record: dict[str, Any]) -> None:
    temporary = path.with_suffix(".tmp")
    with temporary.open("w") as stream:
        stream.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def run_queue(config: dict[str, Any], output: Path) -> None:
    """Wait for all twelve receipts, preserving status and a single command run.

    Args:
        config: Explicit cases, workers, frozen sources, command and wait limits.
        output: New directory reserved for this controller's retained receipts.

    Raises:
        ValueError: Configuration or scientific completion evidence is invalid.
        RuntimeError: An upstream job or the analysis workflow failed.
        TimeoutError: The configured waiting deadline expired.
    """
    output.mkdir(parents=True, exist_ok=False)
    state: dict[str, Any] = {
        "status": "starting",
        "started_utc": _utc(),
        "command": config["command"],
        "cwd": config["cwd"],
    }
    try:
        cases = config["cases"]
        if len(cases) != 12 or len({case["id"] for case in cases}) != 12:
            raise ValueError("Exactly twelve unique matched WT cases are required")
        if not config["worker_statuses"]:
            raise ValueError("Explicit upstream worker status paths required")
        if (
            not math.isfinite(config["poll_seconds"])
            or not math.isfinite(config["deadline_hours"])
            or config["poll_seconds"] < 60
            or config["deadline_hours"] <= 0
        ):
            raise ValueError("Poll >=60 seconds and positive deadline required")
        if not isinstance(config["command"], list) or not config["command"]:
            raise ValueError("Explicit nonempty command argv required")
        if not all(isinstance(arg, str) and arg for arg in config["command"]):
            raise ValueError("Command argv must contain nonempty strings")
        if not Path(config["cwd"]).is_dir():
            raise ValueError("Command cwd must exist")
        _write(output / "config.json", config)
        verify_sources(config)
        reader = ReceiptReader()
        deadline = time.monotonic() + config["deadline_hours"] * 3600
        while True:
            if time.monotonic() >= deadline:
                raise TimeoutError("Matched WT completion deadline expired")
            accepted = inspect_cases(config, reader)
            state.update(
                status="waiting", accepted=accepted, required=12, updated_utc=_utc()
            )
            _write(output / "status.json", state)
            LOGGER.info("Matched WT completion: %d/12", len(accepted))
            if len(accepted) == 12:
                break
            time.sleep(config["poll_seconds"])
        verify_topologies(config)
        verify_sources(config)
        state.update(status="running_analysis", analysis_started_utc=_utc())
        _write(output / "status.json", state)
        result = subprocess.run(config["command"], cwd=config["cwd"], check=False)
        state.update(analysis_ended_utc=_utc(), exit_status=result.returncode)
        if result.returncode != 0:
            raise RuntimeError(f"Analysis command failed with exit {result.returncode}")
        state.update(status="completed", ended_utc=_utc())
        _write(output / "status.json", state)
    except BaseException as error:
        state.update(
            status="failed",
            ended_utc=_utc(),
            error=f"{type(error).__name__}: {error}",
            recovery_hint="Inspect upstream receipts and retained logs; no MD resubmission or reconfiguration was attempted.",
        )
        _write(output / "status.json", state)
        LOGGER.exception("Matched WT analysis queue stopped")
        raise
