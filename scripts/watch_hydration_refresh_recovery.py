#!/usr/bin/env python3
"""Recover failed hydration postprocessing after accepted native segments."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
JOURNAL = ROOT / "output_agent_actions/20261021.jsonl"
GATES = (
    "scc_verified",
    "native_process_passed",
    "finite_temperature_passed",
    "acyl_retained",
)


def read_json(path: Path) -> dict[str, Any]:
    """Read a required JSON object."""
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return data


def sha256(path: Path) -> str:
    """Hash a required file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_hashes(expected: dict[str, str]) -> None:
    """Verify a declared recovery source tree against expected hashes."""
    for filename, digest in expected.items():
        if sha256(Path(filename)) != digest:
            raise ValueError(f"Recovery source changed: {filename}")


def journal_action(phase: str, action: str, result: str = "") -> None:
    """Append a required intent or outcome before continuing recovery work."""
    command = [
        sys.executable,
        str(ROOT / "scripts/record_agent_action.py"),
        "--log",
        str(JOURNAL),
        "--agent",
        "hydration_refresh_recovery",
        "--phase",
        phase,
        "--action",
        action,
        "--permission-mode",
        "escalated",
        "--reason",
        "Repair failed hydration postprocessing after accepted native MD without rerunning native MD",
        "--rollback",
        "Preserve failed receipt and logs; recovery outputs are separate and exclusive",
    ]
    if result:
        command.extend(("--result", result))
    subprocess.run(command, cwd=ROOT, check=True)


def verify_sources(config: dict[str, Any], hook: dict[str, Any]) -> None:
    """Require original source snapshots and hook-time source hashes unchanged."""
    snapshot = BASE / "extension_source_snapshot"
    for relative, expected in config["source_snapshot_sha256"].items():
        if sha256(snapshot / relative) != expected:
            raise ValueError(f"Pinned snapshot changed: {relative}")
    before = hook["source_sha256_before"]
    if before != hook["source_sha256_after"]:
        raise ValueError("Original hook detected source changes during execution")
    for name, expected in before.items():
        if sha256(Path(name)) != expected:
            raise ValueError(f"Original hook source changed: {name}")


def eligible(index: int, config: dict[str, Any]) -> tuple[Path, list[str]] | None:
    """Return the failed hook command only after every native gate passes."""
    name = config["segment_names"][index - 1]
    block = BASE / name
    status_path = block / "status.json"
    if not status_path.exists():
        return None
    status = read_json(status_path)
    native_status = status.get("status")
    if native_status in {"failed", "rejected", "aborted"}:
        raise ValueError(f"{name} native segment ended as {native_status}")
    if native_status != "accepted":
        return None
    process = read_json(block / "process.json")
    acceptance = read_json(block / "scientific_acceptance.json")
    if process.get("exit_code") != 0 or status.get("exit_code") != 0:
        raise ValueError(f"{name} native process did not pass")
    if status.get("predecessor") != (
        "accepted_pilot" if index == 1 else config["segment_names"][index - 2]
    ):
        raise ValueError(f"{name} predecessor mismatch")
    predecessor = "accepted_pilot" if index == 1 else config["segment_names"][index - 2]
    predecessor_restart = (
        BASE / ("equil" if index == 1 else predecessor) / "restart.ncrst"
    )
    if sha256(predecessor_restart) != status.get("start_restart_sha256"):
        raise ValueError(f"{name} restart provenance mismatch")
    if process.get("start_restart_sha256") != status.get("start_restart_sha256"):
        raise ValueError(f"{name} process restart provenance mismatch")
    if not all(acceptance.get(gate) is True for gate in GATES):
        raise ValueError(f"{name} scientific acceptance gate failed")
    boundary = acceptance.get("boundary", {})
    if (
        boundary.get("box_unchanged") is not True
        or not isinstance(boundary.get("frozen_max_displacement_A"), (int, float))
        or boundary["frozen_max_displacement_A"] < 0.0
        or boundary["frozen_max_displacement_A"] > 1e-6
    ):
        raise ValueError(f"{name} boundary acceptance failed")
    hook_path = block / "postprocess.json"
    if not hook_path.exists():
        return None
    hook = read_json(hook_path)
    if hook.get("status") in {"pending", "running"}:
        return None
    verify_sources(config, hook)
    prior_recovery = block / "postprocess_recovery_attempt2.json"
    if prior_recovery.exists():
        prior = read_json(prior_recovery)
        if prior.get("exit_code") == 0:
            return block, []
        raise ValueError(f"{name} isolated recovery attempt failed")
    if hook.get("status") == "completed" and hook.get("exit_code") == 0:
        return block, []
    if hook.get("status") != "failed" or hook.get("exit_code") == 0:
        raise ValueError(f"{name} original postprocess is not terminal failed")
    receipt = block / "postprocess_environment_recovery_watcher.json"
    if receipt.exists():
        raise FileExistsError(f"Recovery already attempted for {block.name}")
    return block, [str(item) for item in hook["argv"]]


def recover(index: int, config: dict[str, Any], recovery: dict[str, Any]) -> bool:
    """Recover one accepted block once; return False while it is not ready."""
    result = eligible(index, config)
    if result is None:
        return False
    block, _original_command = result
    if not _original_command:
        return True
    command = [str(value) for value in recovery["command"]]
    sources = recovery["source_sha256"]
    verify_hashes(sources)
    receipt = block / "postprocess_environment_recovery_watcher.json"
    stdout_path = block / "postprocess_environment_recovery_watcher.stdout.log"
    stderr_path = block / "postprocess_environment_recovery_watcher.stderr.log"
    if any(path.exists() for path in (receipt, stdout_path, stderr_path)):
        raise FileExistsError(f"Recovery output already exists for {block.name}")
    action = f"Run isolated hydration postprocess recovery for {block.name}"
    journal_action("intent", action)
    env = os.environ.copy()
    overrides = {
        key: str(value) for key, value in recovery["environment_overrides"].items()
    }
    env.update(overrides)
    try:
        with stdout_path.open("x") as stdout, stderr_path.open("x") as stderr:
            completed = subprocess.run(
                command, cwd=ROOT, env=env, stdout=stdout, stderr=stderr, check=False
            )
        verify_sources(config, read_json(block / "postprocess.json"))
        verify_hashes(sources)
    except Exception as exc:
        journal_action(
            "outcome", action, f"Recovery failed: {type(exc).__name__}: {exc}"
        )
        raise
    outcome = "succeeded" if completed.returncode == 0 else "failed"
    journal_action(
        "outcome", action, f"Recovery {outcome}; exit code {completed.returncode}"
    )
    payload = {
        "source_receipt": "postprocess.json",
        "command": command,
        "recovery_id": recovery.get("id"),
        "environment_overrides": overrides,
        "source_sha256_before": sources,
        "source_sha256_after": sources,
        "exit_code": completed.returncode,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
        "status": outcome,
    }
    with receipt.open("x") as stream:
        json.dump(payload, stream, indent=2)
        stream.write("\n")
    if completed.returncode:
        raise subprocess.CalledProcessError(completed.returncode, command)
    return True


def watch(
    recovery_config: Path,
    start_index: int = 2,
    poll_seconds: int = 120,
    timeout_seconds: int = 86400,
) -> None:
    """Watch remaining segments and recover terminal hook failures once."""
    config = read_json(BASE / "extension_20ps_config.json")
    recovery = read_json(recovery_config)
    deadline = time.monotonic() + timeout_seconds
    for index in range(start_index, len(config["segment_names"]) + 1):
        while time.monotonic() < deadline:
            if recover(index, config, recovery):
                break
            LOGGER.info(
                "%s not ready; polling in %d seconds",
                config["segment_names"][index - 1],
                poll_seconds,
            )
            time.sleep(min(poll_seconds, max(0, deadline - time.monotonic())))
        else:
            raise TimeoutError(
                f"Timed out waiting for {config['segment_names'][index - 1]}"
            )


def main() -> None:
    """Parse bounded watcher options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recovery-config", type=Path, required=True)
    parser.add_argument("--start-index", type=int, default=2)
    parser.add_argument("--poll-seconds", type=int, default=120)
    parser.add_argument("--timeout-seconds", type=int, default=86400)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    watch(
        args.recovery_config, args.start_index, args.poll_seconds, args.timeout_seconds
    )


if __name__ == "__main__":
    main()
