#!/usr/bin/env python3
"""Execute one root-reviewed hydration replica in a fresh native process."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.hydration_replicas import digest, verify_hashes, write_json  # noqa: E402


def run(plan_path: Path, reviewed_sha256: str) -> dict[str, Any]:
    """Run only after exact plan review; stop at any invalid block or candidate."""
    plan_path = plan_path.resolve(strict=True)
    if digest(plan_path) != reviewed_sha256:
        raise ValueError("Plan SHA-256 does not match the root-reviewed value")
    plan = json.loads(plan_path.read_text())
    if (
        plan.get("status") != "PREPARED_ONLY_NOT_LAUNCHED"
        or plan.get("requires_root_review") is not True
    ):
        raise ValueError("Replica plan is not in reviewable, unlaunched state")
    if plan.get("preflight", {}).get("candidate_found") is not False:
        raise ValueError("Preflight candidate screen is not clear")
    branch = Path(plan["base"]).resolve(strict=True)
    if (branch / "replica_execution.json").exists():
        raise FileExistsError(
            "Replica execution receipt already exists; no in-place retry"
        )
    verify_hashes(plan["pinned_sha256"])
    verify_hashes(
        {
            str(branch / "analysis/source_manifest.json"): plan["preflight"][
                "analysis_source_manifest_sha256"
            ],
            str(branch / "analysis/report.json"): plan["preflight"]["report_sha256"],
            str(branch / "analysis/per_water_summary.json"): plan["preflight"][
                "per_water_summary_sha256"
            ],
        }
    )
    receipt: dict[str, Any] = {
        "status": "running",
        "plan_sha256": reviewed_sha256,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "accepted_new_segments": [],
    }
    write_json(branch / "replica_execution.json", receipt, exclusive=True)
    native = branch / "native_source"
    env = dict(
        os.environ,
        PYTHONPATH=str(native),
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    for name in plan["segment_names"]:
        command = [
            sys.executable,
            "-c",
            "import json; from pathlib import Path; "
            "from nyl12_covalent.deacylation_hydration_extension import run_segment; "
            f"r=run_segment(Path({str(branch)!r}), {name!r}); "
            "print(json.dumps(r, allow_nan=False))",
        ]
        result = subprocess.run(
            command, cwd=branch, env=env, capture_output=True, text=True
        )
        (branch / name / "replica_worker.stdout.log").write_text(result.stdout)
        (branch / name / "replica_worker.stderr.log").write_text(result.stderr)
        if result.returncode != 0:
            raise RuntimeError(
                f"Native segment failed ({name}, exit={result.returncode})"
            )
        status = json.loads((branch / name / "status.json").read_text())
        post = json.loads((branch / name / "postprocess.json").read_text())
        if status.get("status") != "accepted" or status.get("exit_code") != 0:
            raise RuntimeError(f"Native scientific acceptance failed: {name}")
        if post.get("status") != "completed" or post.get("exit_code") != 0:
            raise RuntimeError(f"Required hydration postprocess failed: {name}")
        receipt["accepted_new_segments"].append(name)
        count = int(name.split("_")[1])
        gate = [
            sys.executable,
            "-c",
            "from pathlib import Path; "
            "from nyl12_covalent.hydration_followon_gates import report_candidates; "
            f"print('candidate' if report_candidates(Path({str(branch)!r}), {count}) else 'clear')",
        ]
        screened = subprocess.run(
            gate, cwd=branch, env=env, capture_output=True, text=True, check=True
        )
        if screened.stdout.strip().splitlines()[-1] == "candidate":
            receipt.update(status="HOLD_CANDIDATE_REVIEW", held_after=name)
            break
        write_json(branch / "replica_execution.json", receipt)
    if receipt["status"] == "running":
        receipt["status"] = "completed_25ps_review_required"
    receipt["ended_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(branch / "replica_execution.json", receipt)
    return receipt


def main() -> None:
    """Run the explicitly reviewed branch plan."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--reviewed-sha256", required=True)
    args = parser.parse_args()
    run(args.plan, args.reviewed_sha256)


if __name__ == "__main__":
    main()
