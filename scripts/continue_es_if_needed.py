#!/usr/bin/env python3
"""Continue the unrestrained E-S minimization only if its first run exhausts."""

from __future__ import annotations

import argparse
from datetime import UTC, datetime
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.unrestrained_relay import (  # noqa: E402
    coordinates_valid,
    distances,
    es_input,
    native_converged,
    sha256,
    write_json,
)

LOGGER = logging.getLogger(__name__)


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _wait_for_first_run(process_path: Path, poll_seconds: int) -> dict[str, Any]:
    while True:
        process = json.loads(process_path.read_text(encoding="utf-8"))
        if process.get("ended_utc") is not None:
            return process
        time.sleep(poll_seconds)


def _command(sander: Path, folder: Path, topology: Path, starting: Path) -> list[str]:
    return [
        str(sander),
        "-O",
        "-i",
        str(folder / "minimize.mdin"),
        "-p",
        str(topology),
        "-c",
        str(starting),
        "-o",
        str(folder / "minimized.mdout"),
        "-r",
        str(folder / "optimized.ncrst"),
        "-x",
        str(folder / "optimization.nc"),
        "-inf",
        str(folder / "mdinfo"),
    ]


def continue_if_needed(
    es_dir: Path, continuation_cycles: int, poll_seconds: int
) -> int:
    """Wait for the first run and conditionally launch one continuation.

    Args:
        es_dir: Active unrestrained E-S stage directory.
        continuation_cycles: Maximum additional DL-FIND cycles.
        poll_seconds: Delay between process-record checks.

    Returns:
        Zero when the first run or continuation converges; three when the
        continuation exhausts without convergence.

    Raises:
        FileExistsError: If a continuation was already prepared.
        RuntimeError: If the first or continuation engine process fails.
    """
    continuation = es_dir / "continuation_01"
    status_path = es_dir / "conditional_continuation_status.json"
    if continuation.exists():
        raise FileExistsError(continuation)
    inputs = json.loads((es_dir / "inputs.json").read_text(encoding="utf-8"))
    topology = Path(inputs["topology"])
    process_path = es_dir / "process.json"
    write_json(
        status_path,
        {
            "state": "waiting_for_first_run",
            "continuation_cycles": continuation_cycles,
            "poll_seconds": poll_seconds,
            "first_run_process": str(process_path),
            "updated_utc": _now(),
        },
    )
    first_process = _wait_for_first_run(process_path, poll_seconds)
    if first_process.get("exit_status") != 0:
        raise RuntimeError("First E-S sander process failed")
    if native_converged(es_dir / "minimized.mdout"):
        write_json(
            status_path,
            {
                "state": "first_run_converged",
                "continuation_launched": False,
                "updated_utc": _now(),
            },
        )
        LOGGER.info("First E-S run converged; no continuation launched")
        return 0

    continuation.mkdir()
    starting = continuation / "starting.ncrst"
    shutil.copy2(es_dir / "optimized.ncrst", starting)
    continuation_input = es_input(
        (es_dir / "minimize.mdin").read_text(encoding="utf-8"),
        max_cycles=continuation_cycles,
    )
    (continuation / "minimize.mdin").write_text(continuation_input, encoding="utf-8")
    sander = Path("/data/jpc/code/ambertools26/venv/bin/sander")
    command = _command(sander, continuation, topology, starting)
    source_hashes = {
        "first_run_input": sha256(es_dir / "minimize.mdin"),
        "first_run_output": sha256(es_dir / "minimized.mdout"),
        "starting_restart": sha256(starting),
        "topology": sha256(topology),
        "watcher": sha256(Path(__file__)),
    }
    write_json(
        continuation / "inputs.json",
        {
            "command": command,
            "environment": inputs["environment"],
            "maximum_additional_dlfind_cycles": continuation_cycles,
            "restraints": "none; ntr=0, nmropt=0, all atoms active",
            "source_sha256": source_hashes,
        },
    )
    write_json(
        status_path,
        {
            "state": "continuation_running",
            "continuation": str(continuation),
            "updated_utc": _now(),
        },
    )
    process: dict[str, Any] = {"command": command, "started_utc": _now()}
    write_json(continuation / "process.json", process)
    environment = os.environ.copy()
    environment.update(inputs["environment"])
    with (continuation / "engine.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            cwd=continuation,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    process.update(exit_status=result.returncode, ended_utc=_now())
    write_json(continuation / "process.json", process)
    if result.returncode != 0:
        raise RuntimeError("Conditional E-S continuation sander process failed")

    converged = native_converged(continuation / "minimized.mdout")
    valid = coordinates_valid(topology, continuation / "optimized.ncrst", inputs["box"])
    acceptance = {
        "native_converged": converged,
        "finite_coordinates_and_fixed_box": valid,
        "distances_A": distances(topology, continuation / "optimized.ncrst"),
        "restraints": "none",
    }
    write_json(continuation / "acceptance.json", acceptance)
    state = (
        "continuation_converged" if converged and valid else "continuation_exhausted"
    )
    write_json(
        status_path,
        {
            "state": state,
            "continuation": str(continuation),
            "acceptance": acceptance,
            "updated_utc": _now(),
        },
    )
    LOGGER.info("Conditional continuation finished with state %s", state)
    if not (converged and valid):
        return 3

    endpoint = continuation / "optimized.ncrst"
    path_command = [
        sys.executable,
        str(ROOT / "scripts/run_nyl12_path_inspection.py"),
        "--run",
        str(es_dir.parent / "01_coordinate_drive_to_ti1"),
        "--config",
        str(ROOT / "configs/unrestrained_relay_inspection_20261025.json"),
        "--watch",
        "--poll-seconds",
        "15",
        "--",
        sys.executable,
        str(ROOT / "scripts/run_unrestrained_es_relay.py"),
        "resume-path",
        "--endpoint",
        str(endpoint),
    ]
    write_json(
        status_path,
        {
            "state": "coordinate_driving_running",
            "continuation": str(continuation),
            "accepted_es_endpoint": str(endpoint),
            "accepted_es_endpoint_sha256": sha256(endpoint),
            "path_command": path_command,
            "updated_utc": _now(),
        },
    )
    path_result = subprocess.run(path_command, cwd=ROOT, check=False)
    if path_result.returncode != 0:
        raise RuntimeError("Coordinate-driving workflow failed after continuation")
    write_json(
        status_path,
        {
            "state": "coordinate_driving_complete",
            "continuation": str(continuation),
            "accepted_es_endpoint": str(endpoint),
            "accepted_es_endpoint_sha256": sha256(endpoint),
            "path_command": path_command,
            "path_exit_status": path_result.returncode,
            "updated_utc": _now(),
        },
    )
    return 0


def main() -> int:
    """Parse arguments and run the conditional continuation.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--es-dir", required=True, type=Path)
    parser.add_argument("--continuation-cycles", type=int, default=10000)
    parser.add_argument("--poll-seconds", type=int, default=30)
    args = parser.parse_args()
    if args.continuation_cycles < 1 or args.poll_seconds < 1:
        raise ValueError("Cycle and polling values must be positive")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    return continue_if_needed(
        args.es_dir.resolve(), args.continuation_cycles, args.poll_seconds
    )


if __name__ == "__main__":
    raise SystemExit(main())
