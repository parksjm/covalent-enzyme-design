#!/usr/bin/env python3
"""Run restrained DFTB3/MM minimization and a timed 1 ps NVT heating pilot."""

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
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.restrained_qmmm_md import (  # noqa: E402
    heating_input,
    minimization_input,
    proton_restraint,
)

LOGGER = logging.getLogger(__name__)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def _run_stage(
    sander: Path,
    topology: Path,
    start: Path,
    stage: Path,
    environment: dict[str, str],
) -> tuple[Path, float]:
    command = [
        str(sander),
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
    record: dict[str, Any] = {
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write(stage / "process.json", record)
    started = time.monotonic()
    with (stage / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=stage,
            env={**os.environ, **environment},
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    elapsed = time.monotonic() - started
    record.update(
        exit_status=result.returncode,
        elapsed_seconds=elapsed,
        ended_utc=datetime.now(timezone.utc).isoformat(),
    )
    _write(stage / "process.json", record)
    if result.returncode:
        raise RuntimeError(f"Sander failed in {stage.name}")
    return stage / "restart.ncrst", elapsed


def run(args: argparse.Namespace) -> None:
    """Prepare and execute the two-stage restrained QM/MM pilot."""
    sources = (args.topology, args.start, args.audited_mdin, args.sander)
    if not all(path.is_file() for path in sources):
        raise FileNotFoundError("One or more required source files are missing")
    args.output.mkdir(parents=True, exist_ok=False)
    minimum = args.output / "01_restrained_minimization"
    heating = args.output / "02_nvt_heating_1ps"
    minimum.mkdir()
    heating.mkdir()
    audited = args.audited_mdin.read_text()
    restraint = proton_restraint()
    for stage in (minimum, heating):
        (stage / "proton_coordinate.RST").write_text(restraint)
    (minimum / "run.mdin").write_text(minimization_input(audited))
    (heating / "run.mdin").write_text(heating_input(audited))
    environment = {
        "AMBERHOME": str(args.sander.resolve().parent.parent),
        "LD_LIBRARY_PATH": str(args.sander.resolve().parent.parent / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    metadata = {
        "purpose": "Restrained QM/MM environmental relaxation around transferred-proton endpoint",
        "reaction_coordinate": "d(OG1,HG1)-d(N,HG1), flat bottom 0.60-0.85 A",
        "force_constant_kcal_mol_A2": 50.0,
        "positional_restraints": "none",
        "time_step_fs": 0.5,
        "pilot_heating_ps": 1.0,
        "sources": [
            {"path": str(path.resolve()), "sha256": _sha256(path)} for path in sources
        ],
        "environment": environment,
    }
    _write(args.output / "inputs.json", metadata)
    status: dict[str, Any] = {"state": "running_minimization"}
    _write(args.output / "status.json", status)
    try:
        restart, minimization_seconds = _run_stage(
            args.sander, args.topology, args.start, minimum, environment
        )
        status.update(
            state="running_heating", minimization_seconds=minimization_seconds
        )
        _write(args.output / "status.json", status)
        _, heating_seconds = _run_stage(
            args.sander, args.topology, restart, heating, environment
        )
        status.update(
            state="pilot_completed",
            heating_seconds=heating_seconds,
            estimated_seconds_per_md_step=heating_seconds / 2_000,
        )
    except Exception as error:
        status.update(state="stopped", reason=str(error))
        raise
    finally:
        status["updated_utc"] = datetime.now(timezone.utc).isoformat()
        _write(args.output / "status.json", status)


def main() -> None:
    """Parse arguments and run the restrained QM/MM pilot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--start", required=True, type=Path)
    parser.add_argument("--audited-mdin", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
