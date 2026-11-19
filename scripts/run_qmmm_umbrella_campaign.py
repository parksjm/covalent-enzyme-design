#!/usr/bin/env python3
"""Run direct-PT or asynchronous E-S-to-TI1 QM/MM umbrellas concurrently."""

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


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def window_plan(
    path_run: Path, es_start: Path, coordinate: str
) -> list[dict[str, Any]]:
    """Map the E-S structure and every finalized driven endpoint to umbrellas."""
    if coordinate == "asynchronous_extension":
        status = json.loads((path_run / "extension_status.json").read_text())
        inputs = json.loads((path_run / "extension_inputs.json").read_text())
        labels = status["completed_windows"]
        completed = set(labels)
        schedule = [row for row in inputs["schedule"] if row["label"] in completed]
        if status["state"] != "completed":
            raise ValueError("Asynchronous path extension is not complete")
        if labels != [row["label"] for row in schedule]:
            raise ValueError("Extension schedule and completion order disagree")
        rows = [
            {
                "label": f"window_{int(row['label'].rsplit('_', 1)[1]):02d}",
                "target_A": row["composite_target_A"],
                "start": (
                    path_run
                    / "03_coordinate_drive_extension_to_ti1"
                    / row["label"]
                    / "optimized.ncrst"
                ).resolve(),
            }
            for row in schedule
        ]
        if not rows or not all(row["start"].is_file() for row in rows):
            raise FileNotFoundError("Extension umbrella seed restarts are missing")
        return rows
    status = json.loads((path_run / "status.json").read_text())
    inputs = json.loads((path_run / "inputs.json").read_text())
    labels = status["completed_windows"]
    schedule = inputs["schedule"]
    if labels != [row["label"] for row in schedule]:
        raise ValueError("Coordinate-driving path is not complete and ordered")
    if coordinate == "direct_pt":
        initial_target = inputs["initial_distances_A"]["proton_difference_A"]
        target_key = "target_A"
        drive_directory = "01_coordinate_drive"
    elif coordinate == "asynchronous":
        initial = inputs["initial_geometry_A"]
        initial_target = initial["proton_difference_A"] - initial["OG1_carbonyl_C_A"]
        target_key = "composite_target_A"
        drive_directory = "01_coupled_coordinate_drive"
    else:
        raise ValueError(f"Unsupported umbrella coordinate: {coordinate}")
    rows = [
        {
            "label": "window_00",
            "target_A": initial_target,
            "start": es_start.resolve(),
        }
    ]
    for index, row in enumerate(schedule, start=1):
        rows.append(
            {
                "label": f"window_{index:02d}",
                "target_A": row[target_key],
                "start": (
                    path_run / drive_directory / row["label"] / "optimized.ncrst"
                ).resolve(),
            }
        )
    if not all(row["start"].is_file() for row in rows):
        raise FileNotFoundError("One or more umbrella seed restarts are missing")
    return rows


def stage_plan(protocol: str) -> list[tuple[str, str, int]]:
    """Return the pilot or minimum-30-ps umbrella stage schedule."""
    if protocol == "pilot":
        return [
            ("01_heating_1ps", "heating", 2_000),
            ("02_equilibration_1ps", "equilibration", 2_000),
            ("03_production_2ps", "production", 4_000),
        ]
    if protocol == "production30":
        return [
            ("01_heating_1ps", "heating", 2_000),
            ("02_equilibration_5ps", "equilibration", 10_000),
            *[
                (f"03_production_block_{index:02d}_5ps", "production", 10_000)
                for index in range(1, 7)
            ],
        ]
    raise ValueError(f"Unsupported umbrella protocol: {protocol}")


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
    record: dict[str, Any] = {"command": command, "started_utc": _now()}
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
    protocol: str,
) -> dict[str, Any]:
    folder = output / row["label"]
    folder.mkdir()
    restart = row["start"]
    for name, stage_name, steps in stage_plan(protocol):
        stage = folder / name
        stage.mkdir()
        (stage / "run.mdin").write_text(md_input(audited, stage_name, steps))
        restraint_builder = (
            umbrella_restraint
            if coordinate == "direct_pt"
            else asynchronous_umbrella_restraint
        )
        (stage / "umbrella.RST").write_text(restraint_builder(row["target_A"]))
        restart = _run_stage(stage, sander, topology, restart, environment)
    return {"label": row["label"], "target_A": row["target_A"], "state": "completed"}


def run(args: argparse.Namespace) -> None:
    """Prepare and execute the parallel initial umbrella campaign."""
    sources = (args.topology, args.es_start, args.audited_mdin, args.sander)
    if not all(path.is_file() for path in sources):
        raise FileNotFoundError("Required umbrella source is missing")
    rows = window_plan(args.path_run, args.es_start, args.coordinate)
    if not 1 <= args.workers <= len(rows):
        raise ValueError("workers must be between 1 and the window count")
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
            "purpose": (
                "Initial 300 K DFTB3/MM PMF for direct Thr233-to-N proton transfer"
                if args.coordinate == "direct_pt"
                else "300 K DFTB3/MM PMF extension to bonded TI1 geometry"
                if args.coordinate == "asynchronous_extension"
                else "Initial 300 K DFTB3/MM PMF for asynchronous E-S-to-TI1 coordinate"
            ),
            "reaction_coordinate": (
                "d(OG1,HG1)-d(N,HG1)"
                if args.coordinate == "direct_pt"
                else "d(OG1,HG1)-d(N,HG1)-d(OG1,carbonyl_C)"
            ),
            "temperature_K": 300.0,
            "time_step_fs": 0.5,
            "force_constant_kcal_mol_A2": 50.0,
            "sampling_per_window": {
                "heating_ps": 1.0,
                "equilibration_ps": 1.0 if args.protocol == "pilot" else 5.0,
                "production_ps": 2.0 if args.protocol == "pilot" else 30.0,
                "production_block_ps": (None if args.protocol == "pilot" else 5.0),
            },
            "protocol": args.protocol,
            "workers": args.workers,
            "positional_restraints": "none",
            "windows": [
                {
                    **row,
                    "start": str(row["start"]),
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
                args.protocol,
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
    """Parse campaign inputs and run independent umbrella windows."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path-run", required=True, type=Path)
    parser.add_argument("--es-start", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--audited-mdin", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument(
        "--protocol", choices=("pilot", "production30"), default="pilot"
    )
    parser.add_argument(
        "--coordinate",
        choices=("direct_pt", "asynchronous", "asynchronous_extension"),
        default="direct_pt",
    )
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
