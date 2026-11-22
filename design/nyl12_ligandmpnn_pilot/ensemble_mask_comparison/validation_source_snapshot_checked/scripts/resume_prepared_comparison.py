"""Resume the unchanged prepared comparison after its reviewed preflight repair."""

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

from scripts.run_local_design_validation import run_validation  # noqa: E402


def state(status: str, **details: object) -> None:
    """Record operational state independently of the immutable scientific inputs."""
    path = BASE / "status.json"
    temp = path.with_suffix(".tmp")
    temp.write_text(
        json.dumps(
            {
                "status": status,
                "updated_utc": datetime.now(timezone.utc).isoformat(),
                **details,
            },
            indent=2,
        )
        + "\n"
    )
    os.replace(temp, path)


def main() -> None:
    """Require the reviewed input digest and idle GPUs before restarting preflight."""
    try:
        approval = json.loads(
            (BASE / "prepared/validation_resume_acceptance.json").read_text()
        )
        manifest = BASE / "validation/prepared/candidate_manifest.json"
        if approval["accepted"] is not True or approval["new_prediction_jobs"] != 52:
            raise ValueError("Invalid reviewed resume receipt")
        if (
            hashlib.sha256(manifest.read_bytes()).hexdigest()
            != approval["candidate_manifest_sha256"]
        ):
            raise ValueError("Prepared candidate manifest differs from reviewed input")
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
            state("waiting_for_idle_violet_GPUs", compute_pids=busy.splitlines())
            time.sleep(30)
        else:
            raise TimeoutError("Violet GPUs remained occupied for three hours")
        state(
            "validation_running",
            new_prediction_models=52,
            reused_WT_models=3,
            source="validation_source_snapshot_checked",
            sequence_samples=48,
        )
        result = run_validation(BASE / "validation", BASE.parent / "boltz_validation")
        state(
            "completed",
            validation_completion=result,
            interpretation="E-S structural comparison; catalytic improvement remains untested",
        )
    except Exception as exc:
        state("failed_requires_review", error=f"{type(exc).__name__}: {exc}")
        logging.exception("Prepared comparison stopped; all data are preserved")
        raise


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
