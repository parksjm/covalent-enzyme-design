"""Persist read-only progress snapshots for the active Mistral recovery."""

from __future__ import annotations

import argparse
import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path("/data/jpc/mechanism-agent/benchmarks/nemotron_mistral_nano_20261125")
ATTEMPT = "mistral_immediate_compat_retry01"
LOGGER = logging.getLogger(__name__)


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object from disk."""
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def inspect(root: Path, output: Path) -> dict[str, Any]:
    """Read controller and fresh receipt state without changing job processes."""
    attempt = root / ATTEMPT
    controller = read_json(attempt / "controller_job" / "status.json")
    started = datetime.fromisoformat(controller["started_utc"].replace("Z", "+00:00"))
    cases_doc = read_json(root / "cases.json")
    case_ids = [case["case_id"] for case in cases_doc["cases"]]
    progress_path = root / "results" / "mistral" / "progress.json"
    refs = read_json(progress_path).get("cases", {}) if progress_path.exists() else {}
    successes = 0
    failures = 0
    newest: str | None = None
    for case_id in case_ids:
        ref = refs.get(case_id)
        if not ref:
            continue
        receipt_path = (progress_path.parent / ref["receipt"]).resolve()
        receipt = read_json(receipt_path)
        receipt_started = datetime.fromisoformat(
            receipt["started_utc"].replace("Z", "+00:00")
        )
        if receipt_started < started:
            continue
        ended = receipt.get("ended_utc")
        if ended and (newest is None or ended > newest):
            newest = ended
        if (
            receipt.get("process_status") == "complete"
            and receipt.get("finish_reason") == "stop"
        ):
            successes += 1
        else:
            failures += 1
    phase = (
        "canary/inference"
        if (attempt / "benchmark_resume").exists()
        else "loading_or_warmup"
    )
    exit_status = controller.get("exit_status")
    if exit_status is not None:
        phase = (
            "complete" if exit_status == 0 and successes == len(case_ids) else "failed"
        )
    return {
        "observed_utc": datetime.now(timezone.utc).isoformat(),
        "phase": phase,
        "valid": successes,
        "new_failed": failures,
        "total": len(case_ids),
        "most_recent_completion_utc": newest,
        "controller_exit_status": exit_status,
    }


def persist(output: Path, snapshot: dict[str, Any]) -> None:
    """Write latest snapshot and append an event only when phase/count changes."""
    output.mkdir(parents=True, exist_ok=True)
    latest = output / "latest.json"
    previous = read_json(latest) if latest.exists() else None
    latest.write_text(
        json.dumps(snapshot, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    keys = ("phase", "valid", "new_failed", "total")
    if previous is None or any(previous.get(key) != snapshot[key] for key in keys):
        with (output / "events.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(snapshot, sort_keys=True) + "\n")


def monitor(root: Path, output: Path, poll_seconds: int, max_hours: float) -> int:
    """Poll until controller terminal or timeout; return nonzero for failure."""
    deadline = time.monotonic() + max_hours * 3600
    while True:
        snapshot = inspect(root, output)
        persist(output, snapshot)
        LOGGER.info(
            "phase=%s valid=%d/%d new_failed=%d",
            snapshot["phase"],
            snapshot["valid"],
            snapshot["total"],
            snapshot["new_failed"],
        )
        if snapshot["phase"] in {"complete", "failed"}:
            return 0 if snapshot["phase"] == "complete" else 1
        if time.monotonic() >= deadline:
            snapshot["phase"] = "timeout"
            snapshot["timeout"] = True
            persist(output, snapshot)
            return 2
        time.sleep(poll_seconds)


def main() -> int:
    """Parse CLI options and run the monitor."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--poll", type=int, default=60)
    parser.add_argument("--max-hours", type=float, default=12)
    args = parser.parse_args()
    if args.poll <= 0 or args.max_hours <= 0:
        parser.error("--poll and --max-hours must be positive")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    return monitor(args.root, args.output, args.poll, args.max_hours)


if __name__ == "__main__":
    raise SystemExit(main())
