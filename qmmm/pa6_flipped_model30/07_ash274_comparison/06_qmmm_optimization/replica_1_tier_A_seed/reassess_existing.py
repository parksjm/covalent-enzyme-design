"""Reassess the retained ASH274 endpoint after correcting the log marker."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import runpy
import sys

CASE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "src"))


def sha256(path: Path) -> str:
    """Hash a retained calculation artifact.

    Args:
        path: Existing artifact.

    Returns:
        Hexadecimal SHA-256 digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Write a separate scientific reassessment without rerunning Sander."""
    runner = CASE / "prepare_and_run_ash274_qmmm.py"
    process_path = CASE / "process.json"
    previous_path = CASE / "acceptance.json"
    mdout_path = CASE / "minimized.mdout"
    process = json.loads(process_path.read_text())
    previous = json.loads(previous_path.read_text())
    if process["returncode"] != 0 or process["status"] != "complete":
        raise ValueError("The existing Sander process did not complete")
    if process["output_sha256"]["minimized.mdout"] != sha256(mdout_path):
        raise ValueError("The retained Sander log differs from its process record")
    if previous["acceptance"]["gates"]["normal_engine_completion"]:
        raise ValueError("The previous acceptance has no marker failure to correct")
    output = mdout_path.read_text(errors="replace")
    normal = (
        "Run   done at" in output
        and "Optimisation converged" in output
        and "FATAL" not in output
        and "NaN" not in output
    )
    if not normal:
        raise ValueError("The retained DL-FIND/Sander log lacks completion markers")
    assess = runpy.run_path(str(runner))["assess"]
    acceptance = assess()
    acceptance["gates"]["normal_engine_completion"] = normal
    acceptance["accepted"] = all(acceptance["gates"].values())
    record = {
        "status": "accepted" if acceptance["accepted"] else "rejected",
        "scientific_acceptance": acceptance,
        "reason_for_reassessment": (
            "DL-FIND Sander logs end with 'Run   done at', not 'FINAL RESULTS'"
        ),
        "original_acceptance_sha256": sha256(previous_path),
        "retained_mdout_sha256": sha256(mdout_path),
        "corrected_runner_sha256": sha256(runner),
        "scientific_calculation_rerun": False,
    }
    (CASE / "scientific_acceptance_reassessed.json").write_text(
        json.dumps(record, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
