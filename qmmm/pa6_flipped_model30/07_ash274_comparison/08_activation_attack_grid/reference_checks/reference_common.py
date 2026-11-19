#!/usr/bin/env python3
"""Run the bounded matched-embedding reference queue after baseline and grid completion."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import subprocess
import sys
from typing import Any

HERE = Path(__file__).resolve().parent
GRID = HERE.parent
ROOT = HERE.parents[4]

LOGGER = logging.getLogger(__name__)
TOPOLOGY = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/00_system/ash274_solvated.prmtop"
)
BASELINE = HERE / "es_seed/nonperiodic_full_pc"
BASELINE_JOB = HERE / "job_es_nonperiodic_corrected_route/status.json"
SANDER = Path("/data/jpc/code/ambertools26/venv/bin/sander")
XTB = Path("/data/jpc/code/xtb-6.7.1/xtb")
MAX_WAIT_SECONDS = 48 * 60 * 60


def digest(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Write one readable JSON record."""
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def journal(phase: str, action: str, result: str | None = None) -> None:
    """Append an action record before or after each calculator process."""
    command = [
        sys.executable,
        str(ROOT / "scripts/record_agent_action.py"),
        "--log",
        str(ROOT / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "grid_reference_queue",
        "--phase",
        phase,
        "--permission-mode",
        "escalated",
        "--reason",
        "Run the reviewed bounded matched-embedding reference queue",
        "--rollback",
        "Preserve all completed inputs and logs; use a new attempt directory for retries.",
        "--action",
        action,
    ]
    if result is not None:
        command.extend(("--result", result))
    subprocess.run(command, check=True, cwd=ROOT)
