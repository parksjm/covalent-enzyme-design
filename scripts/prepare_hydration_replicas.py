#!/usr/bin/env python3
"""Prepare independently seeded, review-only hydration branches."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.hydration_replicas import (  # noqa: E402
    BASE_RELATIVE,
    SPECS,
    candidate_check,
    digest,
    prepare_branch,
)

LOGGER = logging.getLogger(__name__)


def prepare_all(base: Path) -> list[dict[str, Any]]:
    """Prepare both branches and screen borrowed history before any MD."""
    base = base.resolve(strict=True)
    results = []
    for spec in SPECS:
        result = prepare_branch(base, spec, sys.executable)
        branch = Path(result["branch"])
        command = [
            sys.executable,
            str(branch / "postprocess_source/scripts/refresh_deacylation_hydration.py"),
            "--base",
            str(branch),
            "--project-root",
            str(branch / "postprocess_source"),
        ]
        env = dict(
            os.environ,
            PYTHONPATH=str(branch / "postprocess_source" / "src"),
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
        )
        subprocess.run(command, cwd=branch, env=env, check=True)
        candidate = candidate_check(
            branch, branch / "native_source", result["checkpoint"], sys.executable
        )
        plan_path = Path(result["plan"])
        plan = json.loads(plan_path.read_text())
        plan["preflight"] = {
            "status": "candidate_found_hold" if candidate else "clear",
            "candidate_found": candidate,
            "analysis_source_manifest_sha256": digest(
                branch / "analysis/source_manifest.json"
            ),
            "report_sha256": digest(branch / "analysis/report.json"),
            "per_water_summary_sha256": digest(
                branch / "analysis/per_water_summary.json"
            ),
        }
        plan["status"] = (
            "HOLD_CANDIDATE_REVIEW" if candidate else "PREPARED_ONLY_NOT_LAUNCHED"
        )
        plan_path.write_text(json.dumps(plan, indent=2, allow_nan=False) + "\n")
        plan_hash = digest(plan_path)
        (branch / "replica_plan.sha256").write_text(plan_hash + "\n")
        result["plan_sha256"] = plan_hash
        result["preflight_candidate"] = candidate
        results.append(result)
    return results


def main() -> None:
    """Run the review-only branch preparation CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / BASE_RELATIVE)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    for item in prepare_all(args.base):
        LOGGER.info(
            "Prepared %s plan=%s sha256=%s candidate=%s",
            item["branch"],
            item["plan"],
            item["plan_sha256"],
            item["preflight_candidate"],
        )


if __name__ == "__main__":
    main()
