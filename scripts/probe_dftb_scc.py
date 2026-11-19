#!/usr/bin/env python3
"""Probe Amber DFTB3 SCC convergence at preserved TS geometries."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import subprocess
from pathlib import Path

LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file.

    Args:
        path: File to hash.

    Returns:
        Hexadecimal digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def probe_input(source: str, scfconv: str, maxiter: int, telec: int) -> str:
    """Build one-gradient input from the original DFTB3/MM settings.

    Args:
        source: Original dimer mdin text.
        scfconv: SCC energy tolerance.
        maxiter: Broyden iterations before reset.
        telec: Electronic temperature in kelvin.

    Returns:
        Single-gradient sander input text.
    """
    qmmm = source.split("&qmmm\n", 1)[1].split("\n/", 1)[0]
    qmmm, count = re.subn(r"scfconv\s*=\s*[^,]+", f"scfconv={scfconv}", qmmm)
    if count != 1:
        raise ValueError("Expected exactly one scfconv setting")
    return (
        "DFTB3/MM one-gradient SCC diagnostic; no dimer or geometry drive\n"
        "&cntrl\n"
        " imin=1, maxcyc=1, ntmin=1, ncyc=0, drms=1.0d-4,\n"
        " ntb=1, igb=0, ntc=1, ntf=1, cut=10.0, ifqnt=1, jfastw=4,\n"
        " ntr=0, ntpr=1, ntwx=0, ntwr=0,\n"
        "/\n&qmmm\n"
        f"{qmmm}\n dftb_maxiter={maxiter}, dftb_telec={telec}, "
        "itrmax=1000, verbosity=1,\n/\n"
    )


def main() -> int:
    """Run controlled SCC tolerance and iteration tests.

    Returns:
        Zero if the diagnostic completed and wrote a summary.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sander", required=True, type=Path)
    parser.add_argument("--only-case", type=str)
    args = parser.parse_args()
    base = args.base.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    topology = base / "inputs/system.prmtop"
    mdin = base / "calculation/dimer.mdin"
    cases = [
        ("start_tight", base / "inputs/ts_guess.ncrst", "1.0d-10", 70, 0),
        ("start_telec_300", base / "inputs/ts_guess.ncrst", "1.0d-8", 70, 300),
        ("checkpoint_tight", base / "calculation/optimized.ncrst", "1.0d-10", 70, 0),
        ("checkpoint_default", base / "calculation/optimized.ncrst", "1.0d-8", 70, 0),
        ("checkpoint_extended", base / "calculation/optimized.ncrst", "1.0d-8", 200, 0),
        (
            "checkpoint_telec_300",
            base / "calculation/optimized.ncrst",
            "1.0d-8",
            70,
            300,
        ),
        (
            "checkpoint_telec_300_tight",
            base / "calculation/optimized.ncrst",
            "1.0d-10",
            70,
            300,
        ),
    ]
    if args.only_case is not None:
        cases = [case for case in cases if case[0] == args.only_case]
        if not cases:
            raise ValueError(f"Unknown diagnostic case: {args.only_case}")
    source = mdin.read_text(encoding="utf-8")
    records: list[dict[str, object]] = []
    for name, coordinates, tolerance, maxiter, telec in cases:
        case_dir = output / name
        case_dir.mkdir()
        input_file = case_dir / "probe.mdin"
        result_file = case_dir / "probe.mdout"
        input_file.write_text(
            probe_input(source, tolerance, maxiter, telec), encoding="utf-8"
        )
        command = [
            str(args.sander),
            "-O",
            "-i",
            str(input_file),
            "-p",
            str(topology),
            "-c",
            str(coordinates),
            "-o",
            str(result_file),
            "-r",
            str(case_dir / "probe.ncrst"),
        ]
        LOGGER.info("Running %s", name)
        with (case_dir / "stdout.log").open("w", encoding="utf-8") as stream:
            result = subprocess.run(
                command,
                cwd=case_dir,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
        mdout = (
            result_file.read_text(encoding="utf-8", errors="replace")
            if result_file.exists()
            else ""
        )
        record = {
            "case": name,
            "command": command,
            "exit_code": result.returncode,
            "scfconv": tolerance,
            "dftb_maxiter": maxiter,
            "dftb_telec_kelvin": telec,
            "scc_failure_count": mdout.count("Convergence could not be achieved"),
            "energy_present": "Energy calculation finished" in mdout,
            "input_sha256": sha256(input_file),
            "coordinate_sha256": sha256(coordinates),
        }
        records.append(record)
        (output / "summary.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )
        LOGGER.info(
            "%s: exit=%s SCC failures=%s",
            name,
            result.returncode,
            record["scc_failure_count"],
        )
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
