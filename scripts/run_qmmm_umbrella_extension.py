#!/usr/bin/env python3
"""Extend a completed QM/MM umbrella pilot in restart-safe production blocks."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import threading
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.qmmm_umbrella import (  # noqa: E402
    asynchronous_umbrella_restraint,
    md_input,
    umbrella_restraint,
)

LOGGER = logging.getLogger(__name__)
TIME_STEP_PS = 0.0005


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def extension_windows(pilot: Path, coordinate: str) -> list[dict[str, Any]]:
    """Return ordered restart inputs from a successfully completed pilot.

    Args:
        pilot: Initial umbrella campaign directory.
        coordinate: Reaction-coordinate identifier.

    Returns:
        Ordered window labels, targets, and production restart paths.

    Raises:
        ValueError: If metadata or coordinate identity is inconsistent.
        FileNotFoundError: If a required restart is absent.
    """
    inputs = json.loads((pilot / "inputs.json").read_text())
    status = json.loads((pilot / "status.json").read_text())
    expected = {
        "direct_pt": "d(OG1,HG1)-d(N,HG1)",
        "asynchronous": "d(OG1,HG1)-d(N,HG1)-d(OG1,carbonyl_C)",
    }
    if coordinate not in expected:
        raise ValueError(f"Unsupported umbrella coordinate: {coordinate}")
    if status.get("state") != "completed" or status.get("failures"):
        raise ValueError("Pilot campaign has not completed successfully")
    if inputs.get("reaction_coordinate") != expected[coordinate]:
        raise ValueError("Pilot reaction-coordinate identity does not match")
    rows = []
    for window in inputs["windows"]:
        label = window["label"]
        if status["windows"].get(label, {}).get("state") != "completed":
            raise ValueError(f"Pilot window is not complete: {label}")
        restart = pilot / label / "03_production_2ps" / "restart.ncrst"
        if not restart.is_file():
            raise FileNotFoundError(restart)
        rows.append({"label": label, "target_A": window["target_A"], "start": restart})
    return rows


def extension_stages(
    equilibration_ps: int, production_ps: int, block_ps: int
) -> list[tuple[str, str, int]]:
    """Build a 5-ps-style equilibration and segmented production schedule."""
    if min(equilibration_ps, production_ps, block_ps) < 1:
        raise ValueError("Stage durations must be positive integer picoseconds")
    if production_ps % block_ps:
        raise ValueError("Production duration must be divisible by block duration")
    stages = [
        (
            f"01_equilibration_{equilibration_ps}ps",
            "equilibration",
            int(equilibration_ps / TIME_STEP_PS),
        )
    ]
    blocks = production_ps // block_ps
    stages.extend(
        (
            f"02_production_block_{index:02d}_{block_ps}ps",
            "production",
            int(block_ps / TIME_STEP_PS),
        )
        for index in range(1, blocks + 1)
    )
    return stages


def _run_stage(
    stage: Path,
    sander: Path,
    topology: Path,
    start: Path,
    environment: dict[str, str],
) -> Path:
    command = [
        str(sander.resolve()),
        "-O",
        "-i",
        str((stage / "run.mdin").resolve()),
        "-p",
        str(topology.resolve()),
        "-c",
        str(start.resolve()),
        "-o",
        str((stage / "run.mdout").resolve()),
        "-r",
        str((stage / "restart.ncrst").resolve()),
        "-x",
        str((stage / "trajectory.nc").resolve()),
        "-inf",
        str((stage / "mdinfo").resolve()),
    ]
    record = {
        "command": command,
        "start_restart_sha256": _sha256(start),
        "started_utc": _now(),
    }
    _write_json(stage / "process.json", record)
    with (stage / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=stage,
            env={**os.environ, **environment},
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(exit_status=result.returncode, ended_utc=_now())
    if result.returncode == 0:
        record["restart_sha256"] = _sha256(stage / "restart.ncrst")
    _write_json(stage / "process.json", record)
    if result.returncode:
        raise RuntimeError(f"Sander failed in {stage}")
    return stage / "restart.ncrst"


def _run_window(
    row: dict[str, Any],
    output: Path,
    topology: Path,
    sander: Path,
    audited: str,
    environment: dict[str, str],
    coordinate: str,
    stages: list[tuple[str, str, int]],
) -> dict[str, Any]:
    folder = output / row["label"]
    folder.mkdir()
    restart = row["start"]
    restraint_builder = (
        umbrella_restraint
        if coordinate == "direct_pt"
        else asynchronous_umbrella_restraint
    )
    for name, stage_name, steps in stages:
        stage = folder / name
        stage.mkdir()
        (stage / "run.mdin").write_text(md_input(audited, stage_name, steps))
        (stage / "umbrella.RST").write_text(restraint_builder(row["target_A"]))
        restart = _run_stage(stage, sander, topology, restart, environment)
    return {"label": row["label"], "target_A": row["target_A"], "state": "completed"}


def run(args: argparse.Namespace) -> None:
    """Validate a pilot and run its equilibration plus production extension."""
    sources = (args.topology, args.audited_mdin, args.sander)
    if not all(path.is_file() for path in sources):
        raise FileNotFoundError("Required umbrella source is missing")
    rows = extension_windows(args.pilot, args.coordinate)
    if not 1 <= args.workers <= len(rows):
        raise ValueError("workers must be between 1 and the window count")
    stages = extension_stages(args.equilibration_ps, args.production_ps, args.block_ps)
    args.output.mkdir(parents=True, exist_ok=False)
    audited = args.audited_mdin.read_text()
    environment = {
        "AMBERHOME": str(args.sander.resolve().parent.parent),
        "LD_LIBRARY_PATH": str(args.sander.resolve().parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    _write_json(
        args.output / "inputs.json",
        {
            "purpose": "Minimum 30-ps DFTB3/MM umbrella production extension",
            "coordinate": args.coordinate,
            "temperature_K": 300.0,
            "time_step_fs": TIME_STEP_PS * 1000.0,
            "force_constant_kcal_mol_A2": 50.0,
            "equilibration_ps": args.equilibration_ps,
            "production_ps": args.production_ps,
            "production_block_ps": args.block_ps,
            "workers": args.workers,
            "positional_restraints": "none",
            "pilot": str(args.pilot.resolve()),
            "pilot_inputs_sha256": _sha256(args.pilot / "inputs.json"),
            "pilot_status_sha256": _sha256(args.pilot / "status.json"),
            "windows": [
                {
                    **row,
                    "start": str(row["start"].resolve()),
                    "start_sha256": _sha256(row["start"]),
                }
                for row in rows
            ],
            "sources": [
                {"path": str(path.resolve()), "sha256": _sha256(path)}
                for path in sources
            ],
        },
    )
    status: dict[str, Any] = {"state": "running", "started_utc": _now(), "windows": {}}
    lock = threading.Lock()
    _write_json(args.output / "status.json", status)
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                _run_window,
                row,
                args.output,
                args.topology,
                args.sander,
                audited,
                environment,
                args.coordinate,
                stages,
            ): row
            for row in rows
        }
        for future in as_completed(futures):
            row = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = {
                    "label": row["label"],
                    "state": "failed",
                    "reason": str(error),
                }
                failures.append(result)
            with lock:
                status["windows"][row["label"]] = result
                status["updated_utc"] = _now()
                _write_json(args.output / "status.json", status)
    status.update(
        state="failed" if failures else "completed",
        failures=failures,
        updated_utc=_now(),
    )
    _write_json(args.output / "status.json", status)
    if failures:
        raise RuntimeError(f"{len(failures)} umbrella windows failed")


def main() -> None:
    """Parse command-line arguments and run the extension."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--audited-mdin", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--equilibration-ps", type=int, default=5)
    parser.add_argument("--production-ps", type=int, default=30)
    parser.add_argument("--block-ps", type=int, default=5)
    parser.add_argument(
        "--coordinate", choices=("direct_pt", "asynchronous"), required=True
    )
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
