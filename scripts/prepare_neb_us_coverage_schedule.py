#!/usr/bin/env python3
"""Freeze a coverage-first controller while retaining the original scientific campaign."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_umbrella_campaign import digest, frozen, write_json  # noqa: E402


def main() -> None:
    """Retain a separately hashed scheduling plan without changing any MD input."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--remote-parent-plan", type=Path, required=True)
    parser.add_argument("--host-python", action="append", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    prepared = args.prepared.resolve(strict=True)
    cfg = frozen(prepared)
    parent_path = prepared.parent / "dispatch/plan.json"
    plan = json.loads(parent_path.read_text())
    remote_parent = json.loads(args.remote_parent_plan.read_text())
    if plan["assignments"] != remote_parent["assignments"]:
        raise ValueError("Host assignments differ between existing controllers")
    hosts = dict(item.split("=", 1) for item in args.host_python)
    if set(hosts) != set(plan["hosts"]):
        raise ValueError("Explicit Python required for each execution host")
    destination = args.output.absolute()
    destination.mkdir(exist_ok=False)
    snapshot = destination / "source_snapshot"
    original = prepared.parent / "source_snapshot"
    for name, expected in cfg["frozen_source_sha256"].items():
        source = Path(name)
        if digest(source) != expected:
            raise ValueError(f"Original scientific source changed: {source}")
        target = snapshot / source.relative_to(original)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    additions = [
        "scripts/reprioritize_neb_umbrella.py",
        "scripts/prepare_neb_us_coverage_schedule.py",
        "src/nyl12_covalent/jobs.py",
        *[
            f"src/nyl12_covalent/neb_umbrella_{part}.py"
            for part in ("coverage", "coverage_worker", "handoff", "segments")
        ],
    ]
    for name in additions:
        target = snapshot / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
    local, remote = plan["local_host"], plan["remote_host"]
    schedule = dict(
        schema_version=1,
        created_utc=datetime.now(timezone.utc).isoformat(),
        authorization="PI: finish active 5ps blocks; prioritize missing windows before extensions; preserve all data",
        hosts=plan["hosts"],
        host_python=hosts,
        local_host=local,
        remote_host=remote,
        windows=plan["windows"],
        assignments=plan["assignments"],
        parent_plans={
            local: str(parent_path),
            remote: str(args.remote_parent_plan.absolute()),
        },
        parent_plan_sha256={
            str(parent_path): digest(parent_path),
            str(args.remote_parent_plan.absolute()): digest(args.remote_parent_plan),
        },
        script=str(snapshot / additions[0]),
        source_sha256={str(p): digest(p) for p in snapshot.rglob("*.py")},
        policy={
            "finish_active_blocks": True,
            "all_63_first_block_gate": True,
            "minimum_production_ps_per_window": 30,
            "production_block_ps": 5,
            "native_inputs_changed": False,
            "native_processes_signaled": False,
            "after_coverage": "Extend least-sampled window in each original worker assignment by one block at a time",
        },
    )
    path = destination / "plan.json"
    write_json(path, schedule)
    path.with_suffix(".sha256").write_text(digest(path) + "\n")
    logging.info("Frozen scheduling-only plan: %s", path)


if __name__ == "__main__":
    main()
