"""Continue the accepted finite design comparison after native sequence sampling."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.design_mask_candidates import prepare_validation  # noqa: E402
from scripts.run_local_design_validation import run_validation  # noqa: E402


def write_state(state: str, **details: object) -> None:
    """Atomically record the current stage without rewriting the frozen protocol."""
    path = BASE / "status.json"
    temp = path.with_suffix(".tmp")
    temp.write_text(
        json.dumps(
            {
                "status": state,
                "updated_utc": datetime.now(timezone.utc).isoformat(),
                **details,
            },
            indent=2,
        )
        + "\n"
    )
    os.replace(temp, path)


def main() -> None:
    """Wait for sampled sequences, prepare matched controls and run local validation."""
    previous = BASE.parent / "boltz_validation"
    generation = BASE / "generation/results/acceptance.json"
    try:
        for _ in range(360):
            if generation.exists():
                receipt = json.loads(generation.read_text())
                required = [
                    "accepted",
                    "native_sampling_completed",
                    "native_sequence_audit_passed",
                    "source_input_hashes_verified",
                    "corrected_export_completed",
                ]
                if (
                    not all(receipt.get(key) is True for key in required)
                    or receipt["samples"] != 48
                ):
                    raise ValueError(f"Sequence acceptance is incomplete: {generation}")
                report = BASE / "generation/results/report.json"
                if (
                    hashlib.sha256(report.read_bytes()).hexdigest()
                    != receipt["report_sha256"]
                ):
                    raise ValueError("Accepted sequence report changed")
                break
            write_state(
                "waiting_for_sequence_generation", generation_job=str(generation)
            )
            time.sleep(30)
        else:
            raise TimeoutError("Sequence generation exceeded three hours")
        # Do not consume GPUs owned by another scientific job at the phase boundary.
        for _ in range(360):
            busy = subprocess.check_output(
                [
                    "nvidia-smi",
                    "--query-compute-apps=pid",
                    "--format=csv,noheader,nounits",
                ],
                text=True,
            ).strip()
            if not busy:
                break
            write_state(
                "waiting_for_idle_violet_GPUs", active_compute_pids=busy.splitlines()
            )
            time.sleep(30)
        else:
            raise TimeoutError("Violet GPUs did not become idle within three hours")
        write_state("preparing_matched_validation")
        report = BASE / "generation/results/report.json"
        prepared = prepare_validation(
            report,
            previous,
            BASE / "validation",
            [2026112201, 2026112202, 2026112203, 2026112204, 2026112205],
        )
        count = len(prepared["jobs"])
        if count > 52:
            raise ValueError(f"Prediction budget exceeded: {count}")
        write_state(
            "validation_running",
            new_prediction_models=count,
            reused_models=len(prepared["reused_models"]),
            sequence_report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
        )
        result = run_validation(BASE / "validation", previous)
        write_state(
            "completed",
            validation_completion=result,
            interpretation="Descriptive E-S structural comparison, not validated catalytic improvement",
        )
    except Exception as exc:
        write_state("failed_requires_review", error=f"{type(exc).__name__}: {exc}")
        logging.exception("Matched comparison stopped; original data are preserved")
        raise


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
