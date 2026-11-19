"""Run a 0 K DFTB3/MM dimer search and stop on the first SCC failure."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import subprocess
import time

LOGGER = logging.getLogger(__name__)
SCC_FAILURE = b"Convergence could not be achieved in this step"


def _digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def run(attempt: Path, sander: Path) -> int:
    """Execute a prepared dimer attempt with explicit process acceptance.

    Args:
        attempt: Prepared attempt directory containing immutable inputs.
        sander: Exact native Amber executable.

    Returns:
        Zero when Sander exits cleanly without an SCC failure; nonzero otherwise.

    Raises:
        FileExistsError: Scientific output already exists in this attempt.
        ValueError: Input provenance or chemistry differs from preparation.
    """
    attempt = attempt.resolve(strict=True)
    sander = sander.resolve(strict=True)
    manifest = json.loads((attempt / "preparation.json").read_text())
    if manifest["status"] != "prepared only; no TS optimization or Hessian run":
        raise ValueError("Attempt is not a prepared-only TS seed")
    if manifest["qm_charge"] != -2 or manifest["qm_atoms"] != 98:
        raise ValueError("QM identity differs from accepted direct NEB")
    inputs = attempt / "inputs"
    calculation = attempt / "calculation"
    for name, expected in manifest["input_sha256"].items():
        path = calculation / name if name == "dimer.mdin" else inputs / name
        if _digest(path) != expected:
            raise ValueError(f"Prepared input changed: {path}")
    mdout = calculation / "dimer.mdout"
    status = calculation / "scientific_status.json"
    if mdout.exists() or status.exists():
        raise FileExistsError("This dimer attempt has already started")
    command = [
        str(sander),
        "-O",
        "-i",
        str(calculation / "dimer.mdin"),
        "-p",
        str(inputs / "system.prmtop"),
        "-c",
        str(inputs / "start.ncrst"),
        "-o",
        str(mdout),
        "-r",
        str(calculation / "optimized.ncrst"),
        "-x",
        str(calculation / "optimization.nc"),
        "-inf",
        str(calculation / "mdinfo"),
    ]
    record: dict[str, object] = {
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "process_exit_status": None,
        "scc_failure": None,
        "native_optimizer_converged": None,
        "first_order_saddle_verified": False,
    }
    status.write_text(json.dumps(record, indent=2) + "\n")
    failure = False
    with (calculation / "stdout.log").open("w") as stream:
        process = subprocess.Popen(
            command, cwd=calculation, stdout=stream, stderr=subprocess.STDOUT
        )
        while process.poll() is None:
            time.sleep(5)
            if mdout.exists() and SCC_FAILURE in mdout.read_bytes():
                failure = True
                LOGGER.error("DFTB3 SCC failure; stopping Sander PID %d", process.pid)
                process.terminate()
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                break
    if mdout.exists() and SCC_FAILURE in mdout.read_bytes():
        failure = True
    output = mdout.read_text(errors="replace") if mdout.exists() else ""
    record.update(
        ended_utc=datetime.now(timezone.utc).isoformat(),
        process_exit_status=process.returncode,
        scc_failure=failure,
        native_optimizer_converged="Optimisation converged" in output,
    )
    status.write_text(json.dumps(record, indent=2) + "\n")
    LOGGER.info("Dimer process ended with %s", record)
    return 0 if process.returncode == 0 and not failure else 1


def main() -> None:
    """Parse the prepared attempt and native Sander executable."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt", type=Path, required=True)
    parser.add_argument("--sander", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(run(args.attempt, args.sander))


if __name__ == "__main__":
    main()
