#!/usr/bin/env python3
"""Run a bounded DL-FIND dimer pilot and stop at the first DFTB SCC failure."""

from __future__ import annotations

import argparse
import json
import logging
import re
import subprocess
import time
from pathlib import Path

LOGGER = logging.getLogger(__name__)
FAILURE_MARKER = b"Convergence could not be achieved in this step"


def pilot_input(source: str, maxcycle: int, mode_path: Path | None) -> str:
    """Change only electronic temperature and pilot cycle limit.

    Args:
        source: Original dimer input.
        maxcycle: Maximum DL-FIND cycles for the bounded pilot.
        mode_path: Absolute displaced-geometry seed for a continuation.

    Returns:
        Pilot input text.
    """
    if "dftb_telec" in source:
        raise ValueError("Original input already specifies dftb_telec")
    source, count = re.subn(
        r"scfconv=1\.0d-10,", "scfconv=1.0d-10, dftb_telec=300,", source
    )
    if count != 1:
        raise ValueError("Expected one original SCC tolerance")
    source, count = re.subn(
        r"tol=1\.d-4, tole=1\.d-5, maxcycle=1000,",
        f"tol=1.d-4, tole=1.d-5, maxcycle={maxcycle},",
        source,
    )
    if count != 1:
        raise ValueError("Expected one DL-FIND cycle limit")
    if mode_path is not None:
        source, count = re.subn(
            r"coords2\(1\)='[^']+'",
            f"coords2(1)='{mode_path.resolve()}'",
            source,
        )
        if count != 1:
            raise ValueError("Expected exactly one dimer-mode input")
    return source


def main() -> int:
    """Run and monitor the pilot, retaining result metadata.

    Returns:
        Zero only if the bounded pilot finishes with no SCC failure.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--maxcycle", type=int, required=True)
    parser.add_argument("--start", type=Path)
    parser.add_argument("--mode", type=Path)
    args = parser.parse_args()
    if (args.start is None) != (args.mode is None):
        raise ValueError("--start and --mode must be provided together")
    if args.maxcycle < 1:
        raise ValueError("--maxcycle must be positive")
    base = args.base.resolve()
    output = args.output.resolve()
    start = args.start.resolve() if args.start else base / "inputs/ts_guess.ncrst"
    mode = args.mode.resolve() if args.mode else None
    if not start.is_file() or (mode is not None and not mode.is_file()):
        raise FileNotFoundError("Dimer start or mode input is missing")
    output.mkdir(parents=True, exist_ok=False)
    original = base / "calculation/dimer.mdin"
    mdin = output / "dimer.mdin"
    mdout = output / "dimer.mdout"
    mdin.write_text(
        pilot_input(original.read_text(encoding="utf-8"), args.maxcycle, mode),
        encoding="utf-8",
    )
    command = [
        str(args.sander),
        "-O",
        "-i",
        str(mdin),
        "-p",
        str(base / "inputs/system.prmtop"),
        "-c",
        str(start),
        "-o",
        str(mdout),
        "-r",
        str(output / "optimized.ncrst"),
        "-x",
        str(output / "optimization.nc"),
        "-inf",
        str(output / "mdinfo"),
    ]
    LOGGER.info("Starting bounded dimer pilot: %s", output)
    failure_seen = False
    with (output / "stdout.log").open("w", encoding="utf-8") as stream:
        proc = subprocess.Popen(
            command, cwd=output, stdout=stream, stderr=subprocess.STDOUT
        )
        while proc.poll() is None:
            time.sleep(2)
            if mdout.exists() and FAILURE_MARKER in mdout.read_bytes():
                failure_seen = True
                LOGGER.error("SCC failure found; terminating sander PID %d", proc.pid)
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                break
    if mdout.exists() and FAILURE_MARKER in mdout.read_bytes():
        failure_seen = True
    record = {
        "command": command,
        "exit_code": proc.returncode,
        "scc_failure": failure_seen,
        "maxcycle": args.maxcycle,
        "dftb_telec_kelvin": 300,
        "scfconv": "1.0d-10",
        "scc_accepted": proc.returncode == 0 and not failure_seen,
    }
    (output / "pilot_status.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info("Pilot ended: %s", record)
    return 0 if record["scc_accepted"] else 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
