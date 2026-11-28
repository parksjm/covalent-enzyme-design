#!/usr/bin/env python3
"""Keep regenerated movie viewers styled for a bounded72-hour period."""

import argparse
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

SOURCE_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("NYL12_VMD_STYLE_REPOSITORY_ROOT", str(SOURCE_ROOT)))
sys.path.insert(0, str(SOURCE_ROOT / "src"))

from nyl12_covalent.vmd_style_watch import inspect_viewer  # noqa: E402

LOGGER = logging.getLogger(__name__)


def run_watch(plan: Path, output: Path, interval: float, hours: float) -> None:
    """Watch logical movie paths, refreshing only stable checksum-ready changes.

    Args:
        plan: Validated explicit movie map with pinned topology hashes.
        output: New job evidence directory for refresh plans/reports/backups.
        interval: Positive poll interval in seconds.
        hours: Positive finite monitoring duration in hours.
    """
    if not 0 < interval <= 60 or not 0 < hours <= 72:
        raise ValueError("Watcher requires interval0-60seconds and duration0-72hours")
    rows = json.loads(plan.read_text())["movies"]
    states: dict[str, dict[str, Any]] = {row["viewer"]: {} for row in rows}
    output.mkdir(parents=True, exist_ok=True)
    end = time.monotonic() + hours * 3600
    cycle = 0
    while time.monotonic() < end:
        cycle += 1
        for row in rows:
            state = states[row["viewer"]]
            if state.get("held"):
                LOGGER.info("HOLD %s: %s", row["viewer"], state["held"])
                continue
            try:
                status = inspect_viewer(row, state)
                LOGGER.info("%s %s", status, row["viewer"])
                if status == "refresh":
                    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
                    event = output / f"viewer_refresh_{stamp}"
                    event.mkdir()
                    explicit = event / "plan.json"
                    explicit.write_text(
                        json.dumps(
                            dict(
                                movies=[
                                    {
                                        **row,
                                        "topology_validation": state.get(
                                            "equivalent_topology"
                                        ),
                                    }
                                ]
                            ),
                            indent=2,
                        )
                        + "\n"
                    )
                    subprocess.run(
                        [
                            sys.executable,
                            str(SOURCE_ROOT / "scripts/refresh_vmd_style.py"),
                            "--plan",
                            str(explicit),
                            "--backup",
                            str(event / "originals"),
                            "--report",
                            str(event / "report.json"),
                            "--apply",
                        ],
                        cwd=ROOT,
                        check=True,
                    )
                    state.clear()
            except (OSError, ValueError, subprocess.CalledProcessError) as exc:
                state["held"] = str(exc)
                LOGGER.error("HOLD %s: %s", row["viewer"], exc)
                subprocess.run(
                    [
                        sys.executable,
                        str(ROOT / "scripts/record_agent_action.py"),
                        "--log",
                        str(ROOT / "output_agent_actions/20261021.jsonl"),
                        "--agent",
                        "vmd_style_watch",
                        "--phase",
                        "outcome",
                        "--permission-mode",
                        "escalated",
                        "--action",
                        f"HOLD display refresh: {row['viewer']}",
                        "--reason",
                        str(exc),
                        "--rollback",
                        "No scientific files touched; inspect exact topology/mapping before restart",
                    ],
                    check=True,
                )
        (output / "watch_status.json").write_text(
            json.dumps(
                dict(
                    cycle=cycle,
                    updated_utc=datetime.now(timezone.utc).isoformat(),
                    states=states,
                ),
                indent=2,
            )
            + "\n"
        )
        time.sleep(min(interval, max(0, end - time.monotonic())))


def main() -> None:
    """Run the bounded view-only monitor with an explicit plan."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval", type=float, default=60)
    parser.add_argument("--hours", type=float, default=72)
    args = parser.parse_args()
    run_watch(args.plan, args.output, args.interval, args.hours)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
