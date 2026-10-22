"""Fail promptly when Sander reports electronically unconverged forces."""

from __future__ import annotations

from pathlib import Path
import subprocess
import time


def run_checked_engine(
    command: list[str], *, cwd: Path, env: dict[str, str], log_path: Path, mdout: Path
) -> dict:
    """Run Sander and stop on its explicit inaccurate-force SCC warning.

    Args:
        command: Exact Sander argument vector.
        cwd: Fresh stage directory.
        env: Explicit subprocess environment.
        log_path: New combined standard output/error log.
        mdout: Engine output to monitor.

    Returns:
        Terminal exit status and whether an SCC failure forced termination.
    """
    warning = "Convergence could not be achieved in this step"
    stopped = False
    with log_path.open("x") as log:
        with subprocess.Popen(
            command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT
        ) as process:
            while process.poll() is None:
                if mdout.exists() and warning in mdout.read_text(errors="replace"):
                    stopped = True
                    process.terminate()
                    try:
                        process.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    break
                time.sleep(0.5)
            code = process.wait()
    if mdout.exists() and warning in mdout.read_text(errors="replace"):
        stopped = True
    return {"exit_status": code, "stopped_on_scf_failure": stopped}
