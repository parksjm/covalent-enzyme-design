"""Nonblocking, block-boundary scheduling for direct-US workers."""

from __future__ import annotations

import json
from contextlib import contextmanager
import fcntl
import logging
import os
from pathlib import Path
import socket
import time
from typing import Any, Iterator

import parmed

from . import direct_us_continuation as continuation
from . import neb_umbrella_campaign as campaign
from .direct_us_continuation_worker import (
    ATTEMPT_FILES,
    _stage_by_block,
    _run_one_stage,
    _worker_event,
)

LOGGER = logging.getLogger(__name__)


@contextmanager
def _window_lock(control: Path, label: str) -> Iterator[bool]:
    """Try to own one window across a single native 5 ps block."""
    lock_dir = control / "locks"
    lock_dir.mkdir(exist_ok=True)
    with (lock_dir / f"{label}.lock").open("a+") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            yield False
            return
        try:
            yield True
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def _window_marker_state(
    plan: dict[str, Any], label: str
) -> tuple[int | None, str | None]:
    """Find first missing block using only durable continuation markers."""
    control = Path(plan["control_dir"])
    blocked_path = control / "receipts" / f"{label}.json"
    if blocked_path.is_file():
        blocked = json.loads(blocked_path.read_text())
        return None, str(blocked.get("reason", "window blocked by prior validation"))
    for block in range(continuation.FIRST_BLOCK, continuation.LAST_BLOCK + 1):
        record = _stage_by_block(plan, label, block)
        stage = Path(record["stage"])
        marker_path = stage / "continuation_status.json"
        if marker_path.is_file():
            marker = json.loads(marker_path.read_text())
            if marker.get("status") == "accepted":
                continue
            return (
                None,
                f"block {block} has nonaccepted attempt state {marker.get('status')}",
            )
        if any((stage / name).exists() for name in ATTEMPT_FILES):
            return (
                None,
                f"block {block} has native attempt artifacts without acceptance",
            )
        return block, None
    return None, None


def _record_blocked(control: Path, label: str, reason: str) -> None:
    """Persist a control-plane block receipt without editing native evidence."""
    path = control / "receipts" / f"{label}.json"
    path.parent.mkdir(exist_ok=True)
    if not path.exists():
        campaign.write_json(
            path, {"status": "blocked", "window": label, "reason": reason}
        )


def _claim_worker(control: Path, worker_id: int, plan_hash: str) -> Path:
    """Create a durable exclusive claim so one worker ID cannot run twice."""
    claim_dir = control / "claims"
    claim_dir.mkdir(exist_ok=True)
    path = claim_dir / f"worker_{worker_id:02d}.json"
    with path.open("x", encoding="utf-8") as stream:
        json.dump(
            {
                "worker_id": worker_id,
                "pid": os.getpid(),
                "host": socket.gethostname(),
                "plan_sha256": plan_hash,
            },
            stream,
            sort_keys=True,
        )
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    return path


def run_worker(
    plan_path: Path, worker_id: int, idle_seconds: float = 2.0
) -> dict[str, Any]:
    """Run dynamic block-boundary work until complete or every window is blocked.

    Args:
        plan_path: Reviewed, immutable continuation plan.
        worker_id: Unique worker number pinned by the plan.
        idle_seconds: Delay before rescanning when other workers hold locks.

    Returns:
        Worker completion or blocked-window summary.
    """
    if type(worker_id) is not int or not 0 <= worker_id < continuation.MAX_WORKERS:
        raise ValueError("worker_id must be an integer from 0 through 39")
    if idle_seconds <= 0:
        raise ValueError("idle_seconds must be positive")
    plan_path = Path(plan_path).resolve(strict=True)
    plan, cfg = continuation.verify_plan(plan_path)
    if worker_id >= plan["workers"]:
        raise ValueError("worker_id exceeds pinned worker count")
    topology = parmed.load_file(cfg["topology"])
    windows = cfg["windows"]
    log_path = plan_path.parent / "workers" / f"worker_{worker_id:02d}.jsonl"
    _claim_worker(plan_path.parent, worker_id, campaign.digest(plan_path))
    log_path.parent.mkdir(exist_ok=True)
    _worker_event(
        log_path,
        "worker_started",
        worker_id=worker_id,
        pid=os.getpid(),
        host=socket.gethostname(),
        plan_sha256=campaign.digest(plan_path),
    )
    blocked: set[str] = set()
    while True:
        made_progress = False
        found_busy = False
        states = {
            row["label"]: _window_marker_state(plan, row["label"]) for row in windows
        }
        for row in windows:
            label = row["label"]
            _, reason = states[label]
            if not reason:
                continue
            with _window_lock(plan_path.parent, label) as acquired:
                if not acquired:
                    found_busy = True
                    continue
                _, confirmed_reason = _window_marker_state(plan, label)
                if not confirmed_reason:
                    found_busy = True
                    continue
                if "nonaccepted attempt state running" in confirmed_reason:
                    confirmed_reason = "abandoned running stage after exclusive window lock became available"
                _record_blocked(plan_path.parent, label, confirmed_reason)
                blocked.add(label)
                _worker_event(
                    log_path,
                    "window_blocked",
                    window=label,
                    reason=confirmed_reason,
                )

        for block in range(continuation.FIRST_BLOCK, continuation.LAST_BLOCK + 1):
            for row in windows:
                label = row["label"]
                candidate, reason = states[label]
                if candidate != block or reason:
                    continue
                with _window_lock(plan_path.parent, label) as acquired:
                    if not acquired:
                        found_busy = True
                        continue
                    current, blocked_reason = _window_marker_state(plan, label)
                    if current != block or blocked_reason:
                        continue
                    try:
                        _run_one_stage(
                            plan_path, plan, cfg, row, topology, block, worker_id
                        )
                    except Exception as error:
                        reason = f"{type(error).__name__}: {error}"
                        _record_blocked(plan_path.parent, label, reason)
                        blocked.add(label)
                        LOGGER.exception(
                            "Continuation %s block %d failed", label, block
                        )
                    made_progress = True
                    break
            if made_progress:
                break
        if made_progress:
            continue

        available = any(
            block is not None and reason is None for block, reason in states.values()
        )
        if not available and not found_busy:
            result = {
                "status": "blocked" if blocked else "completed",
                "blocked_windows": sorted(blocked),
            }
            _worker_event(log_path, "worker_stopped", **result)
            return result
        time.sleep(idle_seconds)
