"""Wait for a validated w033 minimum before invoking a frozen relay driver."""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from nyl12_covalent.relay_path_queue import digest, verify_minimum

LOGGER = logging.getLogger(__name__)


def _event(path: Path, data: dict[str, Any]) -> None:
    with path.open("a") as handle:
        handle.write(
            json.dumps(
                {"utc": datetime.now(timezone.utc).isoformat(), **data}, sort_keys=True
            )
            + "\n"
        )
        handle.flush()
        os.fsync(handle.fileno())


def _journal(phase: str, action: str, files: tuple[Path, ...]) -> None:
    repository = Path(os.environ["NYL12_REPOSITORY"])
    command = [
        sys.executable,
        str(repository / "scripts/record_agent_action.py"),
        "--log",
        str(repository / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "relay_path_watcher",
        "--phase",
        phase,
        "--action",
        action,
        "--permission-mode",
        "escalated",
        "--reason",
        "PI-authorized path conditional on verified w033 minimum",
        "--rollback",
        "Preserve campaign evidence; stop only owned driver process",
    ]
    for path in files:
        command.extend(("--file", str(path)))
    subprocess.run(command, check=True)


def _state(path: Path, data: dict[str, Any]) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def watch(
    verification_root: Path,
    campaign: Path,
    driver: Path,
    poll_seconds: float,
    timeout_hours: float,
) -> int:
    """Wait for the minimum gate and invoke the pinned campaign driver once.

    Args:
        verification_root: Directory containing full-Hessian acceptance.
        campaign: Existing prepared campaign, including config and source snapshot.
        driver: Explicit Python driver inside the campaign source snapshot.
        poll_seconds: Positive delay between gate checks.
        timeout_hours: Positive wall-clock waiting limit.

    Returns:
        The dispatched driver's exit status.

    Raises:
        ValueError: If inputs, acceptance, or pinned artifacts are invalid.
        TimeoutError: If the minimum remains pending until the deadline.
        FileExistsError: If a watcher has already claimed this campaign.
    """
    if poll_seconds <= 0 or timeout_hours <= 0:
        raise ValueError("Polling and timeout must be positive")
    root = verification_root.resolve(strict=True)
    campaign = campaign.resolve(strict=True)
    driver = driver.resolve(strict=True)
    snapshot = (campaign / "source_snapshot").resolve(strict=True)
    if not driver.is_relative_to(snapshot) or not driver.is_file():
        raise ValueError("Driver must be a file in the campaign source snapshot")
    config = campaign / "config.json"
    config_data = json.loads(config.read_text())
    if Path(config_data["verification_root"]).resolve(strict=True) != root:
        raise ValueError("Campaign verification root differs from watcher argument")
    pinned = {"driver_sha256": digest(driver), "config_sha256": digest(config)}
    folder = campaign / "minimum_gate_watcher"
    folder.mkdir(exist_ok=False)
    events, state = folder / "events.jsonl", folder / "wait.json"
    _event(
        events,
        {
            "event": "watch_started",
            "verification_root": str(root),
            "driver": str(driver),
            **pinned,
        },
    )
    deadline = time.monotonic() + timeout_hours * 3600
    try:
        while True:
            if (
                digest(driver) != pinned["driver_sha256"]
                or digest(config) != pinned["config_sha256"]
            ):
                raise ValueError("Pinned campaign driver or config changed")
            receipt = verify_minimum(root)
            if receipt is not None:
                if (
                    digest(driver) != pinned["driver_sha256"]
                    or digest(config) != pinned["config_sha256"]
                ):
                    raise ValueError("Pinned campaign driver or config changed")
                command = [sys.executable, str(driver), "--campaign", str(campaign)]
                data = {
                    "status": "dispatching",
                    "gate": receipt,
                    "command": command,
                    **pinned,
                }
                _state(state, data)
                _event(events, {"event": "dispatch_intent", **data})
                _journal(
                    "intent",
                    f"Dispatch conditional relay driver: {command}",
                    (driver, config, state),
                )
                code = subprocess.call(command)
                _journal(
                    "outcome",
                    f"Conditional relay driver exited {code}: {command}",
                    (driver, config, state),
                )
                _event(events, {"event": "dispatch_outcome", "exit_status": code})
                _state(
                    state,
                    {
                        **data,
                        "status": "completed" if code == 0 else "failed",
                        "exit_status": code,
                    },
                )
                return code
            _state(state, {"status": "waiting", **pinned})
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(
                    "Minimum verification remained pending until timeout"
                )
            time.sleep(min(poll_seconds, remaining))
    except Exception as error:
        _event(events, {"event": "watch_failed", "error": repr(error)})
        _state(state, {"status": "failed", "error": repr(error), **pinned})
        raise


def main() -> None:
    """Parse explicit watcher inputs and run the conditional gate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verification-root", type=Path, required=True)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--poll-seconds", type=float, default=60)
    parser.add_argument("--timeout-hours", type=float, default=24)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    LOGGER.info("Waiting for window033 full-Hessian minimum verification")
    raise SystemExit(
        watch(
            args.verification_root,
            args.campaign,
            args.driver,
            args.poll_seconds,
            args.timeout_hours,
        )
    )


if __name__ == "__main__":
    main()
