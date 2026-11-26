"""Bounded conditional hydration queue; preserves the accepted native chain."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

from .hydration_followon_gates import (
    accepted,
    digest,
    read,
    report_candidates,
    verify_hashes,
)

LOGGER = logging.getLogger(__name__)


def journal(plan: dict[str, Any], phase: str, action: str) -> None:
    """Record intent/outcome, failing closed if append fails."""
    root = Path(plan["root"])
    subprocess.run(
        [
            sys.executable,
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "hydration_followon_queue",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Lead-reviewed bounded hydration continuation",
            "--rollback",
            "Retain all receipts and sources; no automatic retry or deletion",
        ],
        check=True,
    )


def save(base: Path, name: str, payload: dict[str, Any]) -> None:
    """Write a new exclusive finite JSON receipt."""
    with (base / name).open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")


def initial_ready(base: Path) -> bool:
    """Wait for terminal successful current job and complete repaired analysis."""
    job = read(base / "extension_job/status.json")
    if job["ended_utc"] is None:
        return False
    if job["exit_status"] != 0:
        raise ValueError("Original extension process failed")
    status = read(base / "extension_20ps_status.json")
    if status["status"] != "completed" or status["accepted_segments"] != 4:
        raise ValueError("Original extension is not completed/accepted")
    for index in range(1, 5):
        accepted(base, index)
    report = read(base / "analysis/report.json")
    if "extend_04" not in report["frame_count_by_stage"]:
        watcher = read(base / "refresh_recovery_job/status.json")
        if watcher["ended_utc"] is not None and watcher["exit_status"] != 0:
            raise ValueError("Original recovery watcher failed")
        return False
    return True


def hold(plan: dict[str, Any], count: int) -> bool:
    """Save all reported candidate waters for PI review and stop progression."""
    base = Path(plan["base"])
    candidates = report_candidates(base, count)
    if not candidates:
        return False
    save(
        base,
        f"followon_candidate_review_after_{count:02d}.json",
        {
            "status": "PI_review_required",
            "accepted_extension_blocks": count,
            "candidates": candidates,
            "analysis_manifest_sha256": digest(base / "analysis/source_manifest.json"),
            "interpretation": "Plausible geometry only; not a verified catalytic water or reaction.",
        },
    )
    return True


def prepare_runtime(plan: dict[str, Any]) -> None:
    """Create five exclusive native segments from the accepted extend_04 restart."""
    from . import deacylation_hydration as hydration

    base = Path(plan["base"])
    cfg = read(base / "extension_20ps_config.json")
    cfg.update(
        endpoint_restart_path=str(base / "extend_04/restart.ncrst"),
        endpoint_restart_sha256=digest(base / "extend_04/restart.ncrst"),
        segment_names=plan["segment_names"],
        seeds=plan["seeds"],
        postprocess_argv=plan["postprocess_argv"],
        postprocess_environment=plan["postprocess_environment"],
        postprocess_source_sha256=plan["recovery_source_sha256"],
    )
    source = (base / "inputs/source.mdin").read_text()
    for index, name in enumerate(plan["segment_names"]):
        inputs = base / name / "inputs"
        inputs.mkdir(parents=True)
        native = hydration.build_input(
            source, cfg["active_mask"], "equil", plan["seeds"][name]
        )
        (inputs / "segment.mdin").write_text(
            native.replace("hydration only: equil", f"hydration extension: {name}", 1)
        )
        for filename in ("source.mdin", "system.prmtop"):
            (inputs / filename).write_bytes((base / "inputs" / filename).read_bytes())
        if index == 0:
            (inputs / "restart.ncrst").write_bytes(
                (base / "extend_04/restart.ncrst").read_bytes()
            )
        save(
            base / name,
            "status.json",
            {
                "status": "prepared",
                "segment": index + 5,
                "duration_ps": 5.0,
                "predecessor": f"extend_{index + 4:02d}",
            },
        )
    cfg["segment_input_sha256"] = {
        name: {p.name: digest(p) for p in (base / name / "inputs").iterdir()}
        for name in plan["segment_names"]
    }
    snapshot = base / "followon_source_snapshot"
    cfg["source_snapshot_sha256"] = {
        str(p.relative_to(snapshot)): digest(p)
        for p in snapshot.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
    }
    save(base, "followon_runtime_config.json", cfg)


def run(plan_path: Path) -> None:
    """Wait at most 48h; run at most 25ps; stop on candidates or any failed gate."""
    plan = read(plan_path)
    base = Path(plan["base"])
    journal(
        plan,
        "intent",
        f"Conditional bounded queue plan={plan_path}; sha256={digest(plan_path)}",
    )
    try:
        deadline = time.monotonic() + plan["wait_seconds"]
        while True:
            verify_hashes(plan["pinned_sha256"])
            if initial_ready(base):
                break
            if time.monotonic() >= deadline:
                raise TimeoutError("48h current-extension/recovery wait exhausted")
            LOGGER.info("Waiting for accepted 26ps endpoint and repaired full analysis")
            time.sleep(plan["poll_seconds"])
        if hold(plan, 4):
            journal(
                plan,
                "outcome",
                "Candidate found in full 26ps analysis; PI review required; no continuation",
            )
            return
        # A single nice CPU is needed. Preserve evidence before native dispatch.
        resource = subprocess.run(
            ["ps", "-eo", "pid,ppid,ni,pcpu,comm,args"],
            capture_output=True,
            text=True,
            check=True,
        )
        save(
            base,
            "followon_resource_preflight.json",
            {
                "cpu_count": os.cpu_count(),
                "load_average": os.getloadavg(),
                "processes": resource.stdout,
            },
        )
        if os.cpu_count() != 48 or os.getloadavg()[0] > 40:
            raise RuntimeError("CPU resource availability needs review")
        journal(
            plan,
            "intent",
            "Prepare native extend_05..09 only after accepted26ps, pinned sources, no candidate and CPU preflight",
        )
        prepare_runtime(plan)
        journal(
            plan,
            "outcome",
            "Prepared exclusive five blocks from actual accepted extend_04 restart",
        )
        from .deacylation_hydration_extension import run_segment

        for count, name in enumerate(plan["segment_names"], start=5):
            verify_hashes(plan["pinned_sha256"])
            result = run_segment(base, name)
            if (
                result["status"] != "accepted"
                or result["postprocess"]["status"] != "completed"
                or result["postprocess"]["exit_code"] != 0
            ):
                raise ValueError(f"Native acceptance or derived refresh failed: {name}")
            accepted(base, count)
            verify_hashes(plan["pinned_sha256"])
            if hold(plan, count):
                journal(plan, "outcome", f"Candidate after {name}; stop for PI review")
                return
        save(
            base,
            "followon_completion.json",
            {"status": "completed_51ps_review_required", "sampling_converged": False},
        )
        journal(
            plan, "outcome", "Five further blocks accepted: 51ps total; stop for review"
        )
    except Exception as exc:
        journal(
            plan,
            "outcome",
            f"Queue stopped, no automatic retry: {type(exc).__name__}: {exc}",
        )
        raise
