#!/usr/bin/env python3
"""Extend the asynchronous E-S-to-TI1 path until TI1-like geometry is reached."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.coupled_ti1_drive import (  # noqa: E402
    asynchronous_extension_schedule,
    asynchronous_restraint,
    coupled_geometry,
    ti1_like_geometry,
)
from nyl12_covalent.direct_proton_transfer import (  # noqa: E402
    driven_input,
    release_input,
)
from nyl12_covalent.unrestrained_relay import native_converged  # noqa: E402
from scripts.export_direct_pt_progress import export  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def _command(sander: Path, folder: Path, topology: Path, start: Path) -> list[str]:
    folder = folder.resolve()
    return [
        str(sander.resolve()),
        "-O",
        "-i",
        str(folder / "minimize.mdin"),
        "-p",
        str(topology.resolve()),
        "-c",
        str(start.resolve()),
        "-o",
        str(folder / "minimized.mdout"),
        "-r",
        str(folder / "optimized.ncrst"),
        "-x",
        str(folder / "optimization.nc"),
        "-inf",
        str(folder / "mdinfo"),
    ]


def _run_engine(command: list[str], folder: Path, environment: dict[str, str]) -> None:
    record: dict[str, Any] = {"command": command, "started_utc": _now()}
    _write_json(folder / "process.json", record)
    with (folder / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=folder,
            env={**os.environ, **environment},
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(exit_status=result.returncode, ended_utc=_now())
    _write_json(folder / "process.json", record)
    if result.returncode:
        raise RuntimeError(f"Sander failed in {folder}")


def _publish_progress(args: argparse.Namespace) -> Path:
    return export(
        args.path_run,
        args.topology,
        args.cpptraj,
        args.vmd_output,
        "asynchronous",
    )


def run(args: argparse.Namespace) -> None:
    """Continue the completed asynchronous path and test its released endpoint."""
    sources = (
        args.topology,
        args.base_mdin,
        args.sander,
        args.cpptraj,
        args.path_run / "inputs.json",
        args.path_run / "status.json",
    )
    if not all(path.is_file() for path in sources):
        raise FileNotFoundError("Required extension source is missing")
    inputs = json.loads((args.path_run / "inputs.json").read_text())
    status = json.loads((args.path_run / "status.json").read_text())
    schedule = inputs["schedule"]
    if status["state"] != "completed" or len(schedule) < 2:
        raise ValueError("Original asynchronous path must be complete")
    labels = [row["label"] for row in schedule]
    if labels != status["completed_windows"]:
        raise ValueError("Original path schedule and status disagree")
    last_target = float(schedule[-1]["composite_target_A"])
    increment = last_target - float(schedule[-2]["composite_target_A"])
    first_index = int(labels[-1].rsplit("_", 1)[1]) + 1
    rows = asynchronous_extension_schedule(
        last_target, increment, first_index, args.maximum_windows
    )
    windows = args.path_run / "03_coordinate_drive_extension_to_ti1"
    release = args.path_run / "04_unrestrained_extended_ti1_endpoint"
    extension_inputs = args.path_run / "extension_inputs.json"
    extension_status = args.path_run / "extension_status.json"
    for path in (windows, release, extension_inputs, extension_status):
        if path.exists():
            raise FileExistsError(path)
    windows.mkdir()
    release.mkdir()
    start = (
        args.path_run / "01_coupled_coordinate_drive" / labels[-1] / "optimized.ncrst"
    )
    if not start.is_file():
        raise FileNotFoundError(start)
    base = args.base_mdin.read_text()
    environment = {
        "AMBERHOME": str(args.sander.resolve().parent.parent),
        "LD_LIBRARY_PATH": str(args.sander.resolve().parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    _write_json(
        extension_inputs,
        {
            "scientific_intent": "Continue asynchronous driving to a bonded TI1-like endpoint",
            "continuation_start": str(start.resolve()),
            "continuation_start_sha256": _sha256(start),
            "target_increment_A": increment,
            "maximum_windows": args.maximum_windows,
            "early_stop_gate": {
                "OG1_carbonyl_C_max_A": 1.80,
                "N_HG1_max_A": 1.25,
                "OG1_HG1_min_A": 1.30,
                "carbonyl_C_O_min_A": 1.28,
            },
            "schedule": rows,
            "positional_restraints": "none",
            "sources": [
                {"path": str(path.resolve()), "sha256": _sha256(path)}
                for path in sources
            ],
        },
    )
    progress: dict[str, Any] = {
        "state": "running",
        "completed_windows": [],
        "started_utc": _now(),
    }
    _write_json(extension_status, progress)
    restart = start.resolve()
    try:
        for row in rows:
            folder = windows / row["label"]
            folder.mkdir()
            (folder / "minimize.mdin").write_text(driven_input(base))
            (folder / "restraints.RST").write_text(
                asynchronous_restraint(row["composite_target_A"])
            )
            _run_engine(
                _command(args.sander, folder, args.topology, restart),
                folder,
                environment,
            )
            restart = folder / "optimized.ncrst"
            geometry = coupled_geometry(args.topology, restart)
            composite = geometry["proton_difference_A"] - geometry["OG1_carbonyl_C_A"]
            converged = native_converged(folder / "minimized.mdout")
            reached = converged and ti1_like_geometry(geometry)
            _write_json(
                folder / "acceptance.json",
                {
                    "native_converged": converged,
                    "geometry_A": geometry,
                    "composite_coordinate_A": composite,
                    "composite_tracking_error_A": abs(
                        composite - row["composite_target_A"]
                    ),
                    "ti1_like_geometry": reached,
                    "positional_restraints": "none",
                },
            )
            if not converged:
                raise RuntimeError(f"DL-FIND did not converge in {row['label']}")
            progress["completed_windows"].append(row["label"])
            progress.update(latest_geometry_A=geometry, updated_utc=_now())
            _write_json(extension_status, progress)
            bundle = _publish_progress(args)
            LOGGER.info("Published %s", bundle)
            if reached:
                progress["ti1_like_driven_endpoint"] = row["label"]
                break
        else:
            raise RuntimeError("Maximum extension windows reached without TI1 geometry")

        progress.update(state="running_unrestrained_release", updated_utc=_now())
        _write_json(extension_status, progress)
        (release / "minimize.mdin").write_text(release_input(base))
        _run_engine(
            _command(args.sander, release, args.topology, restart),
            release,
            environment,
        )
        final = coupled_geometry(args.topology, release / "optimized.ncrst")
        converged = native_converged(release / "minimized.mdout")
        retained = converged and ti1_like_geometry(final)
        acceptance = {
            "native_converged": converged,
            "ti1_geometry_retained": retained,
            "geometry_A": final,
            "reaction_coordinate_restraints": "none",
            "positional_restraints": "none",
        }
        _write_json(release / "acceptance.json", acceptance)
        progress.update(
            state="completed",
            endpoint_acceptance=acceptance,
            updated_utc=_now(),
        )
        _write_json(extension_status, progress)
        bundle = _publish_progress(args)
        LOGGER.info("Published %s", bundle)
    except Exception as error:
        progress.update(state="stopped", reason=str(error), updated_utc=_now())
        _write_json(extension_status, progress)
        raise


def main() -> None:
    """Parse command-line paths and execute the bounded TI1 extension."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path-run", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--base-mdin", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--vmd-output", required=True, type=Path)
    parser.add_argument("--maximum-windows", type=int, default=8)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
