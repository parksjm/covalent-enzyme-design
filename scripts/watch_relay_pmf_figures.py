#!/usr/bin/env python3
"""Publish standardized uncertainty/overlap figures at cumulative PMF checkpoints."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import time
from zoneinfo import ZoneInfo


def publish_checkpoint(
    checkpoint: Path,
    evidence: Path,
    plots: Path,
    source_root: Path,
    journal_root: Path,
    *,
    replicates: int,
    workers: int,
) -> Path:
    """Compute errors and figures for one immutable completed PMF checkpoint.

    Args:
        checkpoint: Existing analysis checkpoint with JSON report and request.
        evidence: Directory retaining new uncertainty reports and completion records.
        plots: Campaign's existing user-facing plots directory.
        source_root: Frozen source snapshot with scripts and package modules.
        journal_root: Project directory holding the active append-only journal.
        replicates: Explicit five-ps block-bootstrap refit budget.
        workers: Bounded local CPU processes used for the refits.

    Returns:
        Completed publication record, after figures and stable latest copies exist.
    """
    report = checkpoint / "analysis.json"
    request = checkpoint / "analysis_request.json"
    horizon = json.loads(report.read_text())["checkpoint"]["actual_horizon_ps"]
    date = datetime.now(ZoneInfo("America/New_York")).strftime("%Y%m%d")
    stem = f"relay_{horizon:02.0f}ps_standard_{date}"
    completion = evidence / f"checkpoint_{horizon:02.0f}ps_publication.json"
    if completion.exists():
        previous = json.loads(completion.read_text())
        if previous["report_sha256"] != hashlib.sha256(report.read_bytes()).hexdigest():
            raise ValueError(
                "A published checkpoint changed; preserve and review its provenance"
            )
        return completion
    uncertainty = evidence / f"checkpoint_{horizon:02.0f}ps_uncertainty.json"
    journal_command = [
        sys.executable,
        str(journal_root / "scripts/record_agent_action.py"),
        "--log",
        str(journal_root / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "relay_pmf_standard_figures",
        "--action",
        stem,
        "--permission-mode",
        "escalated",
        "--reason",
        "Derived statistical errors and standard figures for a completed checkpoint",
        "--rollback",
        "Preserve all native data; retain derived reports and source snapshot",
        "--file",
        str(report),
        "--file",
        str(request),
    ]
    subprocess.run([*journal_command, "--phase", "intent"], check=True)
    evidence.mkdir(parents=True, exist_ok=True)
    if not uncertainty.exists():
        subprocess.run(
            [
                sys.executable,
                str(source_root / "scripts/analyze_relay_pmf_uncertainty.py"),
                "--report",
                str(report),
                "--request",
                str(request),
                "--output",
                str(uncertainty),
                "--bootstrap-replicates",
                str(replicates),
                "--bootstrap-workers",
                str(workers),
            ],
            check=True,
        )
    subprocess.run(
        [
            sys.executable,
            str(source_root / "scripts/plot_relay_pmf_diagnostics.py"),
            "--report",
            str(report),
            "--uncertainty",
            str(uncertainty),
            "--output-dir",
            str(plots),
            "--stem",
            stem,
            "--interval-method",
            "block-bootstrap",
        ],
        check=True,
    )
    for suffix in ("pmf", "mbar_overlap", "neighbor_overlap"):
        for extension in ("png", "pdf"):
            destination = plots / f"{suffix}_latest.{extension}"
            temporary = destination.with_name(destination.name + ".pending")
            shutil.copy2(plots / f"{stem}_{suffix}.{extension}", temporary)
            temporary.replace(destination)
    completion.write_text(
        json.dumps(
            dict(
                checkpoint=str(checkpoint),
                report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
                uncertainty=str(uncertainty),
                plots=str(plots),
                stem=stem,
                native_jobs_changed=False,
                scientific_acceptance="Statistical diagnostics only; no convergence claim",
            ),
            indent=2,
        )
        + "\n"
    )
    subprocess.run(
        [
            *journal_command,
            "--phase",
            "outcome",
            "--file",
            str(completion),
            "--result",
            "Completed standardized plots and pointwise block-bootstrap intervals",
        ],
        check=True,
    )
    return completion


def main() -> None:
    """Watch existing reports only; never submit, stop, or extend scientific jobs."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("checkpoints", "evidence", "plots", "source-root", "journal-root"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--minimum-horizon", type=int, required=True)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--max-hours", type=float, required=True)
    parser.add_argument("--poll-seconds", type=float, default=60)
    args = parser.parse_args()
    if args.max_hours <= 0 or args.poll_seconds < 30:
        raise ValueError(
            "Positive runtime bound and poll interval at least30s required"
        )
    deadline = time.monotonic() + args.max_hours * 3600
    while time.monotonic() < deadline:
        for checkpoint in sorted(args.checkpoints.glob("checkpoint_*ps")):
            report = checkpoint / "analysis.json"
            request = checkpoint / "analysis_request.json"
            if not report.is_file() or not request.is_file():
                continue
            value = json.loads(report.read_text())
            if value["checkpoint"]["actual_horizon_ps"] < args.minimum_horizon:
                continue
            publish_checkpoint(
                checkpoint,
                args.evidence,
                args.plots,
                args.source_root,
                args.journal_root,
                replicates=args.replicates,
                workers=args.workers,
            )
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
