#!/usr/bin/env python3
"""Watch accepted direct-US continuation prefixes and publish matched PMFs."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import sys
import time

from nyl12_covalent.relay_umbrella_analysis_io import digest

from nyl12_covalent.direct_fixed_watch_helpers import (  # noqa: F401
    HORIZONS_PS,
    POLL_SECONDS,
    _repository_root,
    _record_action,
    _relay_request,
    _continuation_ready,
    _verify_report,
    _journaled_child,
    _plot_metadata,
    _replace_latest,
    _write_status,
)

LOGGER = logging.getLogger(__name__)


def watch(args: argparse.Namespace) -> None:
    """Analyze available direct horizons, then publish only matched relay comparisons."""
    campaign = args.campaign.resolve(strict=True)
    plan = args.continuation_plan.resolve(strict=True)
    relay_evidence = args.relay_evidence.resolve(strict=True)
    source = args.source.resolve(strict=True)
    analyzer = source / "scripts/analyze_direct_fixed_pmf.py"
    plotter = source / "scripts/plot_direct_fixed_pmf.py"
    if not analyzer.is_file() or not plotter.is_file():
        raise FileNotFoundError(
            "Source snapshot must include direct analyzer and plotter scripts"
        )
    root = _repository_root(campaign)
    request = _relay_request(campaign)
    args.evidence.mkdir(parents=True, exist_ok=True)
    args.plots.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + args.max_hours * 3600
    status_path = args.evidence / "direct_fixed_pmf_watch_status.json"
    while time.monotonic() < deadline:
        direct_reports: dict[float, Path] = {}
        for horizon in HORIZONS_PS:
            report_path = args.evidence / f"direct_{horizon:.0f}ps_analysis.json"
            ready, _, plan_sources, augmented = _continuation_ready(
                campaign, plan, horizon
            )
            if report_path.exists():
                _verify_report(report_path, horizon, plan_sources)
                direct_reports[horizon] = report_path
                continue
            if not ready:
                _write_status(
                    status_path,
                    {
                        "status": "waiting_for_direct_checkpoint",
                        "next_horizon_ps": horizon,
                        "updated_utc": datetime.now(timezone.utc).isoformat(),
                        "native_jobs_changed": False,
                    },
                )
                break
            command = [
                sys.executable,
                str(analyzer),
                "--campaign",
                str(campaign),
                "--augmented",
                str(augmented),
                "--relay-request",
                str(request),
                "--output",
                str(report_path),
                "--horizon",
                str(int(horizon)),
                "--replicates",
                str(args.replicates),
                "--workers",
                str(args.workers),
                "--continuation-plan",
                str(plan),
            ]
            _journaled_child(
                root,
                args.evidence,
                f"direct_{horizon:.0f}ps_analysis",
                "PI-authorized fixed-burn-in direct-US continuation analysis",
                command,
                args.evidence / "watcher_logs" / f"direct_{horizon:.0f}ps_analysis.log",
                report_path,
                deadline - time.monotonic(),
                source,
            )
            _verify_report(report_path, horizon, plan_sources)
            direct_reports[horizon] = report_path

        for horizon, direct_report in direct_reports.items():
            relay_report = relay_evidence / f"fixed_{horizon:.0f}ps_analysis.json"
            if not relay_report.is_file():
                _write_status(
                    status_path,
                    {
                        "status": "waiting_for_matching_relay_report",
                        "horizon_ps": horizon,
                        "direct_report": str(direct_report),
                        "relay_report_expected": str(relay_report),
                        "updated_utc": datetime.now(timezone.utc).isoformat(),
                    },
                )
                break
            _verify_report(relay_report, horizon, expected_windows=57)
            plot_output = args.plots
            metadata = _plot_metadata(plot_output, direct_report, relay_report)
            if metadata is None:
                command = [
                    sys.executable,
                    str(plotter),
                    "--direct-report",
                    str(direct_report),
                    "--relay-report",
                    str(relay_report),
                    "--output-dir",
                    str(plot_output),
                ]
                _journaled_child(
                    root,
                    args.evidence,
                    f"direct_{horizon:.0f}ps_matched_plot",
                    "Waited for same-horizon relay report before plotting",
                    command,
                    args.evidence / "watcher_logs" / f"direct_{horizon:.0f}ps_plot.log",
                    plot_output,
                    deadline - time.monotonic(),
                    source,
                )
                metadata = _plot_metadata(plot_output, direct_report, relay_report)
                if metadata is None:
                    raise ValueError("Plot child completed without source metadata")
            publication = args.evidence / f"direct_{horizon:.0f}ps_publication.json"
            if publication.exists():
                previous = json.loads(publication.read_text())
                if previous.get("direct_report_sha256") != digest(
                    direct_report
                ) or previous.get("relay_report_sha256") != digest(relay_report):
                    raise ValueError(
                        f"Existing publication identity changed: {publication}"
                    )
                continue
            aliases = _replace_latest(args.plots, metadata)
            _record_action(
                root,
                "intent",
                f"direct_{horizon:.0f}ps_publish",
                "Publish matched fixed-burn-in direct/relay figure aliases",
                direct_report,
            )
            payload = {
                "horizon_ps": horizon,
                "direct_report": str(direct_report),
                "direct_report_sha256": digest(direct_report),
                "relay_report": str(relay_report),
                "relay_report_sha256": digest(relay_report),
                "plots": metadata["plots"],
                "plot_sha256": {path: digest(Path(path)) for path in metadata["plots"]},
                "latest_aliases": aliases,
                "native_jobs_changed": False,
                "updated_utc": datetime.now(timezone.utc).isoformat(),
            }
            publication.write_text(json.dumps(payload, indent=2) + "\n")
            _record_action(
                root,
                "outcome",
                f"direct_{horizon:.0f}ps_publish",
                "Publish matched fixed-burn-in direct/relay figure aliases",
                publication,
                "Matched horizon published; previous figure targets retained",
            )
        if (args.evidence / "direct_60ps_analysis.json").is_file() and (
            args.evidence / "direct_60ps_publication.json"
        ).is_file():
            LOGGER.info("Direct 60 ps PMF and matched plot are complete")
            return
        time.sleep(min(POLL_SECONDS, max(1.0, deadline - time.monotonic())))
    _write_status(
        status_path,
        {
            "status": "watch_deadline_reached",
            "max_hours": args.max_hours,
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "native_jobs_changed": False,
        },
    )


def main() -> None:
    """Poll accepted continuation checkpoints for at most 72 hours."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--continuation-plan", type=Path, required=True)
    parser.add_argument("--relay-evidence", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--plots", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--max-hours", type=float, default=72.0)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=200)
    args = parser.parse_args()
    if (
        not 0 < args.max_hours <= 72
        or not 1 <= args.workers <= 8
        or args.replicates < 50
    ):
        raise ValueError(
            "Require max-hours in (0,72], workers in 1..8, and positive replicates"
        )
    watch(args)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
