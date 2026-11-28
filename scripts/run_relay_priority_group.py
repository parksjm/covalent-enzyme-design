#!/usr/bin/env python3
"""Run reviewed relay corrections while reserving an exact process group."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import time
from typing import Any

from relay_campaign_resume import identity, matches, resume_all, save

LOGGER = logging.getLogger(__name__)


def verify_files(files: dict[str, str]) -> None:
    """Require every explicitly pinned file to retain its SHA-256 digest."""
    if not files:
        raise ValueError("A nonempty source/input manifest is required")
    for name, expected in files.items():
        path = Path(name)
        if not path.is_absolute():
            raise ValueError("Pinned paths must be absolute")
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Pinned file changed: {name}")


def journal(config: dict[str, Any], phase: str, action: str) -> None:
    """Record process-control intent/outcome before dependent mutations."""
    subprocess.run(
        [
            "/usr/bin/python3",
            config["logger"],
            "--log",
            config["journal"],
            "--agent",
            "relay_priority_group",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "PI prioritized relay PMF repair over direct US and hydration",
            "--rollback",
            "Resume only the recorded process identities; retain data",
        ],
        check=True,
    )


def pause_all(processes: list[dict[str, Any]], boot_id: str) -> list[int]:
    """Pause a fully verified same-owner group using race-safe pidfds."""
    if Path("/proc/sys/kernel/random/boot_id").read_text().strip() != boot_id:
        raise ValueError("Boot changed; refusing process signals")
    descriptors: list[int] = []
    stopped: list[int] = []
    try:
        for row in processes:
            if row["uid"] != os.getuid():
                raise ValueError("Refusing a different owner's process")
            descriptor = os.pidfd_open(row["pid"])
            descriptors.append(descriptor)
            actual = identity(row["pid"])
            if not matches(actual, row) or actual["state"] in {"T", "t", "Z"}:
                raise ValueError("Reserved process identity/state changed")
        for descriptor in descriptors:
            signal.pidfd_send_signal(descriptor, signal.SIGSTOP)
            stopped.append(descriptor)
    except BaseException:
        for descriptor in stopped:
            signal.pidfd_send_signal(descriptor, signal.SIGCONT)
        raise
    finally:
        for descriptor in descriptors:
            os.close(descriptor)
    return [row["pid"] for row in processes]


def preflight(config: dict[str, Any]) -> None:
    """Validate the bounded task manifest before reserving CPU capacity."""
    if not hasattr(os, "pidfd_open") or not hasattr(signal, "pidfd_send_signal"):
        raise RuntimeError("Use the validated system Python pidfd APIs")
    verify_files(config["sha256"])
    tasks = config["tasks"]
    if not 1 <= len(tasks) <= 37 or len({t["label"] for t in tasks}) != len(tasks):
        raise ValueError("Require 1..37 uniquely named reviewed tasks")
    if len(config["processes"]) != 42:
        raise ValueError("Expected the complete 40-direct/2-hydration group")
    for row in tasks:
        if not row["argv"] or not Path(row["argv"][0]).is_absolute():
            raise ValueError("Explicit absolute task executable required")
        if not Path(row["cwd"]).is_dir() or not Path(row["log"]).parent.is_dir():
            raise ValueError("Task working/log directory is missing")
        if Path(row["log"]).exists() or Path(row["receipt"]).exists():
            raise FileExistsError("Priority task already attempted")


def run(config: dict[str, Any]) -> None:
    """Run all corrections concurrently and restore reserved engines afterward.

    A failed task does not terminate other scientific children. The controller
    waits for all launched tasks, then restores the entire reserved group.
    Process completion and scientific acceptance remain separate records.
    """
    preflight(config)
    state = Path(config["status"])
    journal(config, "intent", "Pause all 40 direct and 2 hydration engines together")
    paused = pause_all(config["processes"], config["boot_id"])
    save(state, {"status": "paused", "paused_pids": paused})
    handles: list[Any] = []
    children: list[tuple[dict[str, Any], subprocess.Popen[Any]]] = []
    results: list[dict[str, Any]] = []
    failure: BaseException | None = None
    try:
        journal(config, "outcome", "Paused all 42 verified native engines")
        for row in config["tasks"]:
            verify_files(config["sha256"])
            journal(config, "intent", f"Launch priority relay task {row['label']}")
            handle = Path(row["log"]).open("x")
            handles.append(handle)
            environment = dict(os.environ)
            environment.update(row["environment"])
            child = subprocess.Popen(
                row["argv"],
                cwd=row["cwd"],
                env=environment,
                stdout=handle,
                stderr=subprocess.STDOUT,
            )
            children.append((row, child))
            save(
                state,
                {
                    "status": "running",
                    "paused_pids": paused,
                    "tasks": [{"label": r["label"], "pid": c.pid} for r, c in children],
                },
            )
    except BaseException as error:
        failure = error
    finally:
        for row, child in children:
            code = child.wait()
            receipt_path = Path(row["receipt"])
            try:
                receipt = (
                    json.loads(receipt_path.read_text())
                    if receipt_path.exists()
                    else None
                )
            except (OSError, ValueError) as error:
                receipt = None
                failure = error
                LOGGER.exception("Invalid task receipt: %s", receipt_path)
            accepted = (
                code == 0
                and receipt is not None
                and receipt.get("status") == "completed"
            )
            results.append(
                {
                    "label": row["label"],
                    "exit_status": code,
                    "receipt": str(receipt_path),
                    "accepted": accepted,
                }
            )
            save(
                state,
                {
                    "status": "waiting_for_group",
                    "results": results,
                    "paused_pids": paused,
                },
            )
        for handle in handles:
            handle.close()
        journal(
            config,
            "intent",
            "All launched relay corrections idle; restore reserved group",
        )
        resumed = resume_all(config["processes"], config["boot_id"])
        accepted = (
            failure is None
            and len(results) == len(config["tasks"])
            and all(r["accepted"] for r in results)
        )
        save(
            state,
            {
                "status": "completed" if accepted else "failed",
                "resumed_pids": resumed,
                "tasks": results,
                "error": repr(failure) if failure else None,
                "sampling_converged": False,
                "updated_epoch": time.time(),
            },
        )
        journal(
            config,
            "outcome",
            "Restored all 42 engines; priority task receipts retained",
        )
    if failure is not None:
        raise failure
    if not accepted:
        raise RuntimeError("Priority tasks require review; reserved group was restored")


def main() -> None:
    """Load a caller-pinned plan and run the bounded priority group."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    payload = args.config.read_bytes()
    if hashlib.sha256(payload).hexdigest() != args.sha256:
        raise ValueError("Priority configuration hash changed")
    logging.basicConfig(level=logging.INFO)
    run(json.loads(payload))


if __name__ == "__main__":
    main()
