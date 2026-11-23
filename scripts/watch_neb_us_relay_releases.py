"""Periodically refresh canonical relay movies and stop after a final refresh."""

from __future__ import annotations

import argparse
import fcntl
import json
import logging
import hashlib
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import time

from scripts.export_neb_us_relay_release_vmd import export

ROOT = Path(__file__).resolve().parents[1]
LOGGER = logging.getLogger(__name__)
NAMES = ("ti1_window032", "ti1_window033", "early_window028")


def journal(phase: str, campaign: Path, description: str) -> None:
    """Append a mandatory action record before dependent visualization changes.

    Args:
        phase: Journal intent or outcome.
        campaign: Exact monitored bounded campaign directory.
        description: Nonsecret action details.
    """
    repository = Path(os.environ["NYL12_REPOSITORY"])
    subprocess.run(
        [
            sys.executable,
            str(repository / "scripts/record_agent_action.py"),
            "--log",
            str(repository / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "neb_us_relay_monitor",
            "--phase",
            phase,
            "--action",
            description,
            "--permission-mode",
            "escalated",
            "--reason",
            "Refresh authorized canonical optimizer movies",
            "--rollback",
            "Raw immutable seed inputs and optimizer trajectories retained",
            "--file",
            str(campaign / "prepared.json"),
        ],
        check=True,
    )


def refresh(campaign: Path) -> bool:
    """Publish verified derived movies atomically and report all-jobs completion.

    Args:
        campaign: Bounded three-seed release campaign.

    Returns:
        Whether all three launcher status records have terminal exit statuses.
    """
    journal("intent", campaign, "Refresh three canonical relay optimizer movies")
    with (campaign / ".vmd_export.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            return _refresh_locked(campaign)
        except Exception as error:
            journal(
                "outcome",
                campaign,
                f"Canonical refresh failed: {type(error).__name__}: {error}",
            )
            raise
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def _refresh_locked(campaign: Path) -> bool:
    """Publish all three verified bundles while holding the campaign lock."""
    complete = []
    results = []
    for name in NAMES:
        seed = campaign / name
        status = json.loads((seed / "job/status.json").read_text())
        complete.append(status["exit_status"] is not None)
        temporary = Path(tempfile.mkdtemp(prefix=f"relay_movie_{name}_")) / "bundle"
        export(seed, temporary)
        destination = seed / "vmd_progress_verified"
        destination.mkdir(exist_ok=True)
        for source in sorted(temporary.rglob("*")):
            target = destination / source.relative_to(temporary)
            if source.is_dir():
                target.mkdir(exist_ok=True)
            else:
                pending = target.with_name("." + target.name + ".pending")
                shutil.copy2(source, pending)
                if (
                    hashlib.sha256(source.read_bytes()).digest()
                    != hashlib.sha256(pending.read_bytes()).digest()
                ):
                    raise RuntimeError(f"Movie publication copy mismatch: {source}")
                pending.replace(target)
                source.unlink()
        for directory in sorted(temporary.rglob("*"), reverse=True):
            if directory.is_dir():
                directory.rmdir()
        temporary.rmdir()
        temporary.parent.rmdir()
        manifest = json.loads((destination / "manifest.json").read_text())
        results.append(
            {
                "seed": name,
                "frame_count": manifest["frame_count"],
                "process_terminal": complete[-1],
                "source_snapshot_sha256": manifest["source_snapshot_sha256"],
            }
        )
    checkpoint = {
        "all_jobs_terminal": all(complete),
        "updated_unix": time.time(),
        "seeds": results,
    }
    (campaign / "vmd_refresh_status.json").write_text(
        json.dumps(checkpoint, indent=2) + "\n"
    )
    journal(
        "outcome",
        campaign,
        f"Published three verified canonical movies; all_terminal={all(complete)}",
    )
    LOGGER.info("Updated canonical movies; all_terminal=%s", all(complete))
    return all(complete)


def watch(campaign: Path, interval: float, hours: float) -> None:
    """Refresh immediately, periodically, and after all processes terminate.

    Args:
        campaign: Prepared release campaign.
        interval: Positive interval between refreshes, in seconds.
        hours: Positive maximum monitor lifetime, in hours.
    """
    if interval <= 0 or hours <= 0:
        raise ValueError("Positive interval and hours required")
    deadline = time.monotonic() + hours * 3600
    while True:
        if refresh(campaign):
            return
        if time.monotonic() >= deadline:
            raise RuntimeError(
                "Monitor time limit reached with unfinished scientific jobs"
            )
        time.sleep(interval)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--interval", type=float, default=300)
    parser.add_argument("--hours", type=float, default=168)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    watch(args.campaign.resolve(), args.interval, args.hours)
