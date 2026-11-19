#!/usr/bin/env python3
"""Atomically export accepted relay campaign progress without running calculations."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.relay_progress import (  # noqa: E402
    COMPARISON,
    EXPORT_SOURCE,
    accepted_rows,
    export_movie,
    plot_progress,
    snapshot_frames,
    write_table,
)

LOGGER = logging.getLogger(__name__)
DEFAULT = COMPARISON / "09_water_relay_acylation"


def digest(path: Path) -> str:
    """Return SHA-256 of a source or accepted geometry artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    """Atomically replace a small JSON report within its destination filesystem."""
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def journal(phase: str, action: str) -> None:
    """Journal refresh intent/outcome and fail explicitly on logging failures."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_campaign_progress",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "tool",
            "--reason",
            "Accepted-only stage-aware campaign plots/movie; no energy evaluation",
            "--rollback",
            "Generated immutable progress updates can be regenerated; preserve source snapshots",
        ],
        check=True,
    )


def finished(campaign: Path) -> bool:
    """Return explicit campaign terminal status, independent of chemical acceptance."""
    path = campaign / "campaign_status.json"
    if not path.exists():
        return False
    status = json.loads(path.read_text())
    return status.get("finished") is True or status.get("status") in {
        "finished",
        "completed",
        "failed",
        "stopped",
        "blocked",
    }


def refresh(campaign: Path, output: Path) -> bool:
    """Publish one complete immutable update and atomically switch its pointer.

    Args:
        campaign: Directory containing case.json and native stage subdirectories.
        output: Progress directory; prior updates and original artifacts persist.

    Returns:
        Whether campaign_status.json declares a terminal process status.
    """
    case_path = campaign / "case.json"
    case = json.loads(case_path.read_text())
    rows = accepted_rows(campaign)
    terminal = finished(campaign)
    output.mkdir(parents=True, exist_ok=True)
    stability = output / "snapshot_stability.json"
    previous = json.loads(stability.read_text()) if stability.is_file() else {}
    observed, stable = {}, []
    for path, row in snapshot_frames(campaign, rows):
        stat = path.stat()
        stamp = [stat.st_size, stat.st_mtime_ns]
        observed[str(path)] = stamp
        if terminal or previous.get(str(path)) == stamp:
            stable.append((path, row))
    atomic_json(stability, observed)
    sources = [
        Path(__file__).resolve(),
        ROOT / "src/nyl12_covalent/relay_progress.py",
        ROOT / "src/nyl12_covalent/relay_progress_vmd.py",
        EXPORT_SOURCE / "export_water_path.py",
        EXPORT_SOURCE / "water_path_vmd_support.py",
    ]
    ledgers = sorted(
        {campaign / row["stage"] / "accepted_iterations.jsonl" for row in rows}
    )
    identity = {
        "case_sha256": digest(case_path),
        "source_sha256": {str(p): digest(p) for p in sources},
        "accepted_jsonl_sha256": {str(p): digest(p) for p in ledgers},
        "stable_frame_sha256": {str(p): digest(p) for p, _ in stable},
        "terminal": terminal,
    }
    report_path = output / "progress.json"
    if (
        report_path.is_file()
        and json.loads(report_path.read_text()).get("identity") == identity
    ):
        return terminal
    number = len(list(output.glob("update_[0-9]*"))) + 1
    name = f"update_{number:06d}"
    pending = output / f".building_{number:06d}"
    pending.mkdir()  # Interrupted builds remain evidence; never overwrite one.
    source_folder = pending / "source"
    source_folder.mkdir()
    for path in sources:
        shutil.copy2(path, source_folder / path.name)
    atomic_json(source_folder / "manifest.json", identity)
    atomic_json(pending / "accepted_iterations.json", rows)
    if rows:
        write_table(pending / "accepted_iterations.csv", rows)
        plot_progress(pending, rows)
        export_movie(pending, stable, case)
    report = {
        "remote_evidence_sha256": {
            "accepted_iterations.json": digest(pending / "accepted_iterations.json")
        },
        "identity": identity,
        "accepted_records": len(rows),
        "terminal": terminal,
        "latest_update": name,
        "stable_stage_snapshots": len(stable),
        "frame_definition": "Accepted stage initial/every25/final only; identical geometry references share movie frames",
        "interpretation": "Primary-water screen negative; independent accepted-ES endpoint test; no relay-ready, PMF or validated-barrier claim",
        "energy_reference": "First accepted physical energy; bias shown absolute",
    }
    atomic_json(pending / "progress.json", report)
    payload = sorted(
        p
        for p in pending.rglob("*")
        if p.is_file()
        and not p.is_symlink()
        and "frame_links" not in p.relative_to(pending).parts
        and p.name != "accepted_iterations.json"
    )
    (pending / "MANIFEST.sha256").write_text(
        "".join(f"{digest(path)}  {path.relative_to(pending)}\n" for path in payload)
    )
    pending.rename(output / name)
    pointer = output / f".latest_{os.getpid()}"
    pointer.symlink_to(name, target_is_directory=True)
    pointer.replace(output / "latest")
    atomic_json(report_path, report)
    LOGGER.info(
        "Published %s: accepted=%d snapshots=%d terminal=%s",
        name,
        len(rows),
        len(stable),
        terminal,
    )
    return terminal


def main() -> None:
    """Refresh once or monitor every30s until explicit campaign completion."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=DEFAULT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=float, default=30.0)
    parser.add_argument("--hours", type=float, default=48.0)
    args = parser.parse_args()
    if not (0 < args.interval <= 60) or not 0 < args.hours < float("inf"):
        parser.error("interval must be finite in (0,60] seconds; hours finite positive")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    campaign = args.campaign.resolve()
    output = args.output.resolve() if args.output else campaign / "progress"
    deadline = time.monotonic() + 3600 * args.hours
    while True:
        journal("intent", f"Refresh accepted campaign progress: {campaign} -> {output}")
        try:
            terminal = refresh(campaign, output)
        except Exception as error:
            journal(
                "outcome", f"Progress refresh failed: {type(error).__name__}: {error}"
            )
            raise
        journal("outcome", f"Progress refresh completed; terminal={terminal}")
        if terminal or not args.watch:
            return
        if time.monotonic() >= deadline:
            raise TimeoutError("Progress monitor reached bounded deadline")
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
