#!/usr/bin/env python3
"""Audit an authorized read-only RAG run without submitting or changing jobs."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Callable

LOGGER = logging.getLogger(__name__)
SCIENCE_ROOT = Path("/data/jpc/covalent-enzyme-design")
AGENT_ROOT = Path("/data/jpc/mechanism-agent")


def record(phase: str, detail: str) -> None:
    """Append a checked action record to both established journals.

    Args:
        phase: Journal phase accepted by the project logger.
        detail: Nonsecret description of this monitoring action.
    """
    for journal in (
        SCIENCE_ROOT / "output_agent_actions/20261021.jsonl",
        AGENT_ROOT / "output_agent_actions/actions.jsonl",
    ):
        subprocess.run(
            [
                sys.executable,
                str(SCIENCE_ROOT / "scripts/record_agent_action.py"),
                "--log", str(journal), "--agent", "rag_overnight_audit",
                "--phase", phase, "--action", "observe_readonly_pilot",
                "--permission-mode", "escalated", "--reason", detail,
                "--rollback", "Read-only observations; preserve jobs and raw evidence",
            ],
            check=True,
        )


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Atomically publish a monitoring report without touching native records.

    Args:
        path: Monitor-owned destination.
        data: JSON-safe report content.
    """
    temporary = path.with_suffix(path.suffix + ".pending")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


def notify(event: str, run: Path, message: str) -> None:
    """Use the existing trusted recipient/configuration for one factual alert.

    Args:
        event: Existing notification event type.
        run: Audit directory identified in the notification.
        message: Counts and limitations derived from retained artifacts.
    """
    subprocess.run(
        [
            str(AGENT_ROOT / ".venv/bin/mechanism-agent"), "notify",
            "--notifications-config", str(AGENT_ROOT / "deployment/notifications.toml"),
            "--event", event, "--run", str(run), "--message", message,
        ],
        check=True,
    )


def needs_followup(native: dict[str, Any], summary: dict[str, Any]) -> bool:
    """Identify an incomplete/error run needing a factual correction alert.

    Args:
        native: Native controller record, kept separate from scientific acceptance.
        summary: Actual response/citation counts from the audited pilot.

    Returns:
        Whether the ordinary completion email is insufficient or inaccurate.
    """
    return (
        native["exit_status"] != 0
        or summary["response_count"] != summary["expected_calls"]
        or summary["citation_valid_count"] != summary["expected_calls"]
        or bool(summary["issues"])
    )


def monitor(
    pilot: Path,
    controller: Path,
    output: Path,
    summarize: Callable[[Path, Path], dict[str, Any]],
    format_message: Callable[[dict[str, Any]], str],
    *,
    hours: float = 12.0,
    interval: float = 60.0,
    sleep: Callable[[float], None] = time.sleep,
    monotonic: Callable[[], float] = time.monotonic,
) -> None:
    """Wait for terminal receipts, then retain a truthful morning audit.

    This monitor never edits inputs, restarts a model, or submits a scientific
    calculation. Success/failure notifications from the existing controller are
    supplemented at most once when actual counts or validation need attention.

    Args:
        pilot: Authorized pilot root.
        controller: Exact current controller status file.
        output: New monitor-owned state directory.
        summarize: Read-only reporting helper.
        format_message: Formatter for actual artifact counts and limitations.
        hours: Bounded monitoring lifetime.
        interval: Poll interval in seconds, at least one second.
        sleep: Injectable wait for focused tests.
        monotonic: Injectable deadline clock.
    """
    if hours <= 0 or interval < 1:
        raise ValueError("Positive duration and poll interval >= 1 second required")
    output.mkdir(parents=True, exist_ok=False)
    record("intent", "PI-authorized overnight RAG audit; no job dispatch or relaunch")
    deadline = monotonic() + hours * 3600
    while True:
        record("intent", "Read current pilot receipts and actual saved response counts")
        summary = summarize(pilot, controller)
        native = json.loads(controller.read_text())
        observed = {
            "observed_utc": datetime.now(timezone.utc).isoformat(),
            "native_controller": str(controller),
            "native_exit_status": native["exit_status"],
            "summary": summary,
            "scientific_submission_enabled": False,
        }
        write_json(output / "latest.json", observed)
        with (output / "events.jsonl").open("a") as stream:
            stream.write(json.dumps(observed) + "\n")
        record("outcome", "Monitor snapshot saved; no scientific or model jobs changed")
        terminal = native["ended_utc"] is not None and native["exit_status"] is not None
        if terminal or monotonic() >= deadline:
            timed_out = not terminal
            observed["monitor_deadline_reached"] = timed_out
            message = format_message(summary)
            if timed_out:
                message = "Audit deadline reached; controller was left unchanged. " + message
            write_json(output / "final_report.json", observed)
            (output / "MORNING_REPORT.md").write_text(
                "# Read-only RAG overnight audit\n\n"
                + message + "\n\n"
                + "Factual acceptance requires human/lead review. No scientific jobs "
                "were submitted, stopped, or modified by this monitor.\n"
            )
            if timed_out or needs_followup(native, summary):
                event = "question" if timed_out else "error"
                record("intent", "Send one factual follow-up with actual pilot counts")
                notify(event, output, message)
                record("outcome", "Factual follow-up notification command completed")
            record("outcome", "Overnight audit retained; human factual review pending")
            return
        sleep(interval)


def main() -> None:
    """Run the fixed-scope monitor with its frozen sibling reporting helper."""
    from reporting import format_notification, summarize_pilot

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--controller-status", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hours", type=float, default=12.0)
    parser.add_argument("--interval", type=float, default=60.0)
    args = parser.parse_args()
    try:
        monitor(
            args.pilot, args.controller_status, args.output,
            summarize_pilot, format_notification,
            hours=args.hours, interval=args.interval,
        )
    except Exception as error:
        LOGGER.exception("Overnight RAG audit failed; underlying jobs left unchanged")
        record("intent", "Report an audit failure without changing the pilot")
        notify(
            "error", args.output,
            f"The read-only workflow audit failed: {type(error).__name__}: {error}. "
            "Pilot/scientific jobs were left unchanged; consult retained audit logs. "
            "No claim about completed answer count or factual acceptance is made.",
        )
        record("outcome", "Audit failure notification command completed")
        raise


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
