#!/usr/bin/env python3
"""Publish continuous-window PMFs and matched-CI comparisons at five-ps checkpoints."""

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


def publish(
    checkpoint: Path, evidence: Path, plots: Path, source: Path, journal_root: Path
) -> None:
    """Analyze and publish one complete checkpoint without accessing native jobs.

    Args:
        checkpoint: Immutable source report/request directory.
        evidence: Directory for separate continuous analyses and publication records.
        plots: Existing user-facing campaign plots directory.
        source: Immutable source snapshot directory.
        journal_root: Project containing the active append-only journal.
    """
    parent = checkpoint / "analysis.json"
    request = checkpoint / "analysis_request.json"
    horizon = json.loads(parent.read_text())["checkpoint"]["actual_horizon_ps"]
    identity = hashlib.sha256(parent.read_bytes()).hexdigest()
    publication = evidence / f"continuous_{horizon:02.0f}ps_publication.json"
    if publication.exists():
        if json.loads(publication.read_text())["parent_report_sha256"] != identity:
            raise ValueError("Published source checkpoint changed; review provenance")
        return
    journal = [
        sys.executable,
        str(journal_root / "scripts/record_agent_action.py"),
        "--log",
        str(journal_root / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "continuous_relay_pmf_watcher",
        "--action",
        f"continuous_{horizon:g}ps",
        "--permission-mode",
        "escalated",
        "--reason",
        "PI-authorized continuous-window reanalysis",
        "--rollback",
        "Retain derived evidence; original PMFs and native jobs unchanged",
        "--file",
        str(parent),
        "--file",
        str(request),
    ]
    subprocess.run([*journal, "--phase", "intent"], check=True)
    evidence.mkdir(parents=True, exist_ok=True)
    prefix = evidence / f"continuous_{horizon:02.0f}ps"
    report = prefix.with_name(prefix.name + "_analysis.json")
    uncertainty = prefix.with_name(prefix.name + "_uncertainty.json")
    if not report.exists() and not uncertainty.exists():
        subprocess.run(
            [
                sys.executable,
                str(source / "scripts/analyze_relay_continuous_pmf.py"),
                "--request",
                str(request),
                "--parent-report",
                str(parent),
                "--output-stem",
                str(prefix),
                "--replicates",
                "200",
                "--workers",
                "2",
            ],
            check=True,
        )
    if not report.is_file() or not uncertainty.is_file():
        raise ValueError("Partial continuous analysis requires reviewed recovery")
    baseline = evidence / f"blockwise_{horizon:02.0f}ps_uncertainty.json"
    if not baseline.exists():
        subprocess.run(
            [
                sys.executable,
                str(source / "scripts/analyze_relay_pmf_uncertainty.py"),
                "--request",
                str(request),
                "--report",
                str(parent),
                "--output",
                str(baseline),
            ],
            check=True,
        )
    date = datetime.now(ZoneInfo("America/New_York")).strftime("%Y%m%d")
    stem = f"relay_{horizon:02.0f}ps_continuous_{date}"
    subprocess.run(
        [
            sys.executable,
            str(source / "scripts/plot_relay_pmf_diagnostics.py"),
            "--report",
            str(report),
            "--uncertainty",
            str(uncertainty),
            "--output-dir",
            str(plots),
            "--stem",
            stem,
            "--interval-method",
            "analytic",
        ],
        check=True,
    )
    comparison = plots / f"comparison_{horizon:02.0f}ps_continuous_vs_blockwise_{date}"
    subprocess.run(
        [
            sys.executable,
            str(source / "scripts/compare_relay_pmf_selection.py"),
            "--baseline",
            str(baseline),
            "--continuous",
            str(uncertainty),
            "--output-stem",
            str(comparison),
            "--interval-method",
            "analytic",
        ],
        check=True,
    )
    for suffix in ("pmf", "mbar_overlap", "neighbor_overlap"):
        for extension in ("png", "pdf"):
            destination = plots / f"{suffix}_latest.{extension}"
            temporary = destination.with_name(destination.name + ".pending")
            shutil.copy2(plots / f"{stem}_{suffix}.{extension}", temporary)
            temporary.replace(destination)
    publication.write_text(
        json.dumps(
            dict(
                parent_report_sha256=identity,
                report=str(report),
                uncertainty=str(uncertainty),
                comparison=str(comparison),
                plots_stem=stem,
                selection_policy="continuous_window",
                interval_method="Matched PyMBAR asymptotic pointwise95% intervals",
                native_jobs_changed=False,
                sampling_converged=False,
            ),
            indent=2,
        )
        + "\n"
    )
    subprocess.run(
        [
            *journal,
            "--phase",
            "outcome",
            "--file",
            str(publication),
            "--result",
            "Published continuous-window PMF and matched statistical-CI comparison",
        ],
        check=True,
    )


def main() -> None:
    """Watch existing completed checkpoints within a bounded runtime."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("checkpoints", "evidence", "plots", "source", "journal-root"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--minimum-horizon", required=True, type=float)
    parser.add_argument("--max-hours", required=True, type=float)
    args = parser.parse_args()
    if args.max_hours <= 0 or args.minimum_horizon <= 0:
        raise ValueError("Positive horizon and runtime bound required")
    deadline = time.monotonic() + args.max_hours * 3600
    while time.monotonic() < deadline:
        for checkpoint in sorted(args.checkpoints.glob("checkpoint_*ps")):
            parent = checkpoint / "analysis.json"
            if (
                not parent.is_file()
                or not (checkpoint / "analysis_request.json").is_file()
            ):
                continue
            if (
                json.loads(parent.read_text())["checkpoint"]["actual_horizon_ps"]
                >= args.minimum_horizon
            ):
                publish(
                    checkpoint,
                    args.evidence,
                    args.plots,
                    args.source,
                    args.journal_root,
                )
        time.sleep(60)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
