#!/usr/bin/env python3
"""Publish the fixed-burn-in relay PMF at successive complete checkpoints."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import time

from nyl12_covalent.relay_umbrella_analysis_io import digest


def publish(
    checkpoint: Path,
    evidence: Path,
    plots: Path,
    source: Path,
    root: Path,
    workers: int = 2,
) -> None:
    """Refit a complete checkpoint once and atomically update its inspection alias.

    Args:
        checkpoint: Immutable source-request and parent-report directory.
        evidence: Existing logical campaign analysis directory.
        plots: Existing user-facing plot directory.
        source: Frozen analysis/plot source tree.
        root: Repository containing the mandatory action journal.
        workers: Bounded number of bootstrap workers on the analysis host.
    """
    parent = checkpoint / "analysis.json"
    horizon = float(json.loads(parent.read_text())["checkpoint"]["actual_horizon_ps"])
    publication = evidence / f"fixed_{horizon:02.0f}ps_publication.json"
    if publication.exists():
        if json.loads(publication.read_text())["parent_sha256"] != digest(parent):
            raise ValueError("Previously published parent changed")
        return
    journal = [
        sys.executable,
        str(root / "scripts/record_agent_action.py"),
        "--log",
        str(root / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "fixed_burnin_pmf_watcher",
        "--action",
        f"fixed_{horizon:g}ps",
        "--permission-mode",
        "escalated",
        "--reason",
        "PI-authorized fixed-burn-in PMF",
        "--rollback",
        "Retain derived evidence; native jobs unchanged",
    ]
    subprocess.run([*journal, "--phase", "intent", "--file", str(parent)], check=True)
    report = evidence / f"fixed_{horizon:02.0f}ps_analysis.json"
    if not report.exists():
        subprocess.run(
            [
                sys.executable,
                str(source / "scripts/analyze_relay_fixed_pmf.py"),
                "--checkpoint",
                str(checkpoint),
                "--output",
                str(report),
                "--replicates",
                "200",
                "--workers",
                str(workers),
            ],
            check=True,
        )
    if json.loads(report.read_text())["source_sha256"][str(parent.resolve())] != digest(
        parent
    ):
        raise ValueError("Retained analysis refers to another source report")
    from importlib.util import module_from_spec, spec_from_file_location

    specification = spec_from_file_location(
        "fixed_plot", source / "scripts/plot_relay_fixed_pmf.py"
    )
    if specification is None or specification.loader is None:
        raise ImportError("Frozen plot module cannot be loaded")
    plotting = module_from_spec(specification)
    specification.loader.exec_module(plotting)
    metadata = plotting.plot_report(report, plots)
    for filename in metadata["plots"]:
        if "_pmf." not in filename:
            continue
        path = Path(filename)
        alias = plots / f"pmf_latest{path.suffix}"
        pending = alias.with_name(alias.name + ".fixed_pending")
        shutil.copy2(path, pending)
        pending.replace(alias)
    record = dict(
        parent_sha256=digest(parent),
        report=str(report),
        report_sha256=digest(report),
        horizon_ps=horizon,
        plots=metadata["plots"],
        native_jobs_changed=False,
        protocol="fixed5ps_production_burnin_all_frames_contiguous_bootstrap",
    )
    publication.write_text(json.dumps(record, indent=2) + "\n")
    subprocess.run(
        [
            *journal,
            "--phase",
            "outcome",
            "--file",
            str(publication),
            "--result",
            "Published fixed-burn-in PMF, block intervals and sensitivity",
        ],
        check=True,
    )


def main() -> None:
    """Watch complete source checkpoints for a bounded72hours without MD access."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("checkpoints", "evidence", "plots", "source", "journal-root"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--minimum-horizon", type=float, required=True)
    parser.add_argument("--max-hours", type=float, required=True)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.max_hours <= 0 or args.minimum_horizon < 20 or not 1 <= args.workers <= 8:
        raise ValueError("Positive runtime and >=20ps horizon required")
    deadline = time.monotonic() + args.max_hours * 3600
    args.evidence.mkdir(parents=True, exist_ok=True)
    while time.monotonic() < deadline:
        for checkpoint in sorted(args.checkpoints.glob("checkpoint_*ps")):
            parent = checkpoint / "analysis.json"
            if parent.is_file():
                horizon = json.loads(parent.read_text())["checkpoint"][
                    "actual_horizon_ps"
                ]
                if horizon >= args.minimum_horizon:
                    publish(
                        checkpoint,
                        args.evidence,
                        args.plots,
                        args.source,
                        args.journal_root,
                        args.workers,
                    )
        time.sleep(45)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
