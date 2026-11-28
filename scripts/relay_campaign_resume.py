"""Restore an explicitly paused campaign after its priority checks settle."""

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

LOGGER = logging.getLogger(__name__)


def identity(pid: int) -> dict[str, Any]:
    """Read a Linux process identity, including its non-reusable start tick."""
    path = Path("/proc") / str(pid)
    fields = (path / "stat").read_text().split(") ", 1)[1].split()
    return {
        "pid": pid,
        "start_ticks": int(fields[19]),
        "state": fields[0],
        "uid": path.stat().st_uid,
        "argv": (path / "cmdline").read_bytes().replace(b"\0", b" ").decode(),
        "cwd": str((path / "cwd").resolve(strict=True)),
    }


def matches(actual: dict[str, Any], expected: dict[str, Any]) -> bool:
    """Require the recorded PID, start tick, owner, command and directory."""
    return all(
        actual[key] == expected[key]
        for key in ("pid", "start_ticks", "uid", "argv", "cwd")
    )


def live(expected: dict[str, Any]) -> bool:
    """Check whether the exact recorded wrapper is still alive."""
    try:
        actual = identity(expected["pid"])
    except (FileNotFoundError, ProcessLookupError):
        return False
    return matches(actual, expected) and actual["state"] != "Z"


def target_busy(root: Path) -> bool:
    """Detect retained child calculations before releasing reserved capacity."""
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            if process.stat().st_uid != os.getuid():
                continue
            directory = (process / "cwd").resolve(strict=True)
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if directory.is_relative_to(root):
            return True
    return False


def save(path: Path, data: dict[str, Any]) -> None:
    """Publish an atomic operational receipt."""
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


def journal(config: dict[str, Any], phase: str, action: str) -> None:
    """Require an audit record before process-control mutations."""
    subprocess.run(
        [
            "/usr/bin/python3",
            config["logger"],
            "--log",
            config["journal"],
            "--agent",
            "relay_campaign_resume",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "PI-approved relay priority; restore the complete direct campaign together",
            "--rollback",
            "Retain pause identities and every scientific artifact; signal only exact PIDs",
        ],
        check=True,
    )


def resume_all(processes: list[dict[str, Any]], boot_id: str) -> list[int]:
    """Validate all reserved processes before resuming the complete group."""
    if Path("/proc/sys/kernel/random/boot_id").read_text().strip() != boot_id:
        raise ValueError("Boot changed; refusing process signals")
    descriptors: list[int] = []
    try:
        for row in processes:
            if row["uid"] != os.getuid():
                raise ValueError("Refusing to signal a different owner")
            descriptors.append(os.pidfd_open(row["pid"]))
            if not matches(identity(row["pid"]), row):
                raise ValueError("Reserved process identity changed")
        for descriptor in descriptors:
            signal.pidfd_send_signal(descriptor, signal.SIGCONT)
    finally:
        for descriptor in descriptors:
            os.close(descriptor)
    return [row["pid"] for row in processes]


def run(config: dict[str, Any]) -> None:
    """Wait for all four checks, then restore all paused engines together."""
    if not hasattr(os, "pidfd_open") or not hasattr(signal, "pidfd_send_signal"):
        raise RuntimeError(
            "This controller requires the validated system Python pidfd APIs"
        )
    state_path = Path(config["status"])
    journal(
        config, "intent", "Monitor four bounded checks before campaign-wide restoration"
    )
    while True:
        targets = []
        for row in config["targets"]:
            busy = live(row["wrapper"]) or target_busy(Path(row["root"]))
            path = Path(row["root"]) / "window_status.json"
            native = json.loads(path.read_text()) if path.exists() else None
            targets.append({"label": row["label"], "busy": busy, "native": native})
        save(
            state_path,
            {"status": "waiting", "targets": targets, "updated_epoch": time.time()},
        )
        if not any(row["busy"] for row in targets):
            break
        time.sleep(config["poll_seconds"])
    journal(
        config,
        "intent",
        "All priority native checks idle; resume the 40 direct and 2 hydration engines",
    )
    pids = resume_all(config["processes"], config["boot_id"])
    result = {
        "status": "resumed",
        "resumed_pids": pids,
        "targets": targets,
        "updated_epoch": time.time(),
        "scientific_convergence_claimed": False,
        "known_wrapper_issue": config["known_wrapper_issue"],
    }
    save(state_path, result)
    journal(
        config,
        "outcome",
        "Restored all 42 verified native engines together; consult individual scientific acceptance records",
    )
    LOGGER.info("Resumed %d exact native engines together", len(pids))


def main() -> None:
    """Load the pinned control manifest and execute its bounded monitor."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    payload = args.config.read_bytes()
    if hashlib.sha256(payload).hexdigest() != args.sha256:
        raise ValueError("Controller configuration hash changed")
    run(json.loads(payload))


if __name__ == "__main__":
    main()
