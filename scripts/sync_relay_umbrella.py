#!/usr/bin/env python3
"""Bounded compact-only remote relay artifacts and process observations."""

from __future__ import annotations

import argparse
import fcntl
import json
import logging
from pathlib import Path
import shlex
import subprocess
import time


def journal(root: Path, phase: str, action: str) -> None:
    """Append an action before and after each remote observation/compact transfer."""
    subprocess.run(
        [
            "python3",
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_compact_sync",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Authorized bounded native campaign status and compact analysis transfer",
            "--rollback",
            "No remote scientific mutation; keep local copies and logs",
        ],
        check=True,
    )


def sync_once(args: argparse.Namespace) -> None:
    """Refresh status and copy explicitly allowed small artifacts under a marker."""
    marker = args.campaign / ".compact_sync_in_progress"
    marker.write_text(json.dumps(dict(started_epoch=time.time())) + "\n")
    ssh = [
        "ssh",
        "-i",
        str(args.identity),
        "-o",
        "IdentitiesOnly=yes",
        "-o",
        "BatchMode=yes",
    ]
    journal(
        args.journal_root,
        "intent",
        "Observe native relay processes and rsync compact outputs from " + args.remote,
    )
    command = [
        args.remote_python,
        str(args.inspector),
        "--campaign",
        str(args.campaign),
        "--execution",
        str(args.execution),
    ]
    subprocess.run([*ssh, args.remote, shlex.join(command)], check=True)
    include = [
        "*/",
        "*.json",
        "*.jsonl",
        "*.csv",
        "*.npz",
        "*.png",
        "*.sha256",
        "*.mdin",
        "*.RST",
        "umbrella_coordinate.dat",
    ]
    rsync = [
        "rsync",
        "-a",
        "--prune-empty-dirs",
        "--exclude",
        "failed_analysis_*/***",
        "-e",
        shlex.join(ssh),
    ]
    for pattern in include:
        rsync.extend(["--include", pattern])
    rsync.extend(
        [
            "--exclude",
            "*",
            args.remote + ":" + str(args.campaign) + "/",
            str(args.campaign) + "/",
        ]
    )
    subprocess.run(rsync, check=True)
    if args.benchmark_source is not None:
        destination = args.execution / "resource_benchmark"
        destination.mkdir(exist_ok=True)
        subprocess.run(
            [
                "rsync",
                "-a",
                "--prune-empty-dirs",
                "-e",
                shlex.join(ssh),
                "--include",
                "*/",
                "--include",
                "*.json",
                "--exclude",
                "*",
                args.remote + ":" + str(args.benchmark_source) + "/",
                str(destination) + "/",
            ],
            check=True,
        )
    journal(
        args.journal_root,
        "outcome",
        "Native relay status and compact outputs synced successfully from "
        + args.remote,
    )
    marker.unlink()


def main() -> None:
    """Run one exclusive bounded synchronizer, keeping the marker on failures."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("campaign", "execution", "inspector", "identity", "journal-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--benchmark-source", type=Path)
    parser.add_argument("--remote", required=True)
    parser.add_argument("--remote-python", required=True)
    parser.add_argument("--interval", type=float, required=True)
    parser.add_argument("--duration-hours", type=float, required=True)
    args = parser.parse_args()
    if args.interval < 30 or args.duration_hours <= 0:
        raise ValueError("Require bounded positive duration and >=30second interval")
    with (args.execution / "compact_sync.lock").open("a") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        deadline = time.monotonic() + args.duration_hours * 3600
        while time.monotonic() < deadline:
            sync_once(args)
            time.sleep(args.interval)


if __name__ == "__main__":
    main()
