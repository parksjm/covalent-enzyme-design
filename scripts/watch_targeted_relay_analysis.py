#!/usr/bin/env python3
"""Wait at most 24 hours for four accepted target tests, then analyze only."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import subprocess
import sys
import time
from pathlib import Path


def read(path: Path) -> dict:
    """Read a required JSON receipt.

    Args:
        path: Receipt path.

    Returns:
        Parsed receipt.
    """
    return json.loads(path.read_text())


def stage_ready(stage: Path, steps: int) -> bool:
    """Require accepted native completion regardless of outer wrapper exit.

    Args:
        stage: Native stage directory.
        steps: Expected integrated steps.

    Returns:
        Whether the immutable native receipts allow derived analysis.
    """
    required = ("acceptance.json", "process.json", "analysis_outcome.json")
    if not all((stage / n).is_file() for n in required):
        return False
    a, p = read(stage / "acceptance.json"), read(stage / "process.json")
    return (
        p.get("exit_status") == 0
        and bool(p.get("ended_utc"))
        and a.get("accepted") is True
        and a.get("completed_steps") == steps
        and a.get("finite_energy_temperature") is True
        and a.get("every_integration_step_SCC_converged") is True
        and bool(a.get("checkpoint_sha256"))
    )


def workers_alive(target: Path) -> bool:
    """Inspect local process identities without sending signals.

    Args:
        target: Absolute target campaign directory.

    Returns:
        Whether any native engine or target wrapper references this campaign.
    """
    for receipt in target.glob("*/wrapper_state/wrapper_status.json"):
        identity = read(receipt)
        proc = Path("/proc") / str(identity["pid"])
        try:
            fields = (proc / "stat").read_text().rsplit(")", 1)[1].split()
            if int(fields[19]) == identity["start_ticks"] and fields[0] != "Z":
                return True
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            command = (
                (entry / "cmdline")
                .read_bytes()
                .replace(b"\x00", b" ")
                .decode(errors="replace")
            )
            cwd = (entry / "cwd").resolve(strict=True)
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if ("sander" in command and target in cwd.parents) or (
            "run_relay_targeted.py" in command and str(target) in command
        ):
            return True
    return False


def journal(root: Path, phase: str, action: str) -> None:
    """Append mandatory execution evidence before dependent mutations.

    Args:
        root: Repository journal root.
        phase: Intent or outcome.
        action: Nonsecret description.
    """
    subprocess.run(
        [
            sys.executable,
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "targeted_relay_analysis",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "PI-authorized post-run diagnostic analysis",
            "--rollback",
            "Retain derived evidence; frozen native inputs untouched",
        ],
        check=True,
    )


def run(
    prepared: Path, output: Path, root: Path, max_hours: float, poll: float
) -> None:
    """Collect accepted stages and emit success or explicit partial/failure packet.

    Args:
        prepared: Frozen target preparation JSON with adjacent checksum.
        output: Existing targeted analysis directory.
        root: Action journal repository.
        max_hours: Bounded wait, up to 24 hours.
        poll: Receipt polling interval in seconds.
    """
    import parmed
    from nyl12_covalent.relay_umbrella_runtime import ensure_diagnostics
    from relay_targeted_report import make_report

    if not 0 < max_hours <= 24 or poll < 10 or not output.is_dir():
        raise ValueError("Require existing output, 0<runtime<=24h and poll>=10s")
    expected = prepared.with_suffix(".sha256").read_text().strip()
    if hashlib.sha256(prepared.read_bytes()).hexdigest() != expected:
        raise ValueError("Frozen preparation checksum changed")
    cfg = read(prepared)
    for name, checksum in cfg["frozen_files_sha256"].items():
        if (
            hashlib.sha256((prepared.parent / name).read_bytes()).hexdigest()
            != checksum
        ):
            raise ValueError(f"Frozen target input changed: {name}")
    topology = None
    completed = set()
    started = time.monotonic()
    journal(
        root, "intent", "Bounded post-run watcher; derived-only accepted-stage analysis"
    )
    try:
        while time.monotonic() - started < max_hours * 3600:
            failed = []
            all_done = True
            for row in cfg["windows"]:
                status_path = prepared.parent / row["label"] / "window_status.json"
                status = read(status_path) if status_path.exists() else {}
                for spec in row["stages"]:
                    stage = (prepared.parent / spec["input"]).parent
                    if stage in completed:
                        continue
                    if not stage_ready(stage, spec["steps"]):
                        all_done = False
                        continue
                    journal(root, "intent", f"Ensure derived diagnostics {stage}")
                    if topology is None:
                        topology = parmed.load_file(cfg["topology"])
                    if not ensure_diagnostics(cfg, row, stage, topology):
                        raise RuntimeError(f"Derived diagnostics failed: {stage}")
                    completed.add(stage)
                    journal(root, "outcome", f"Derived diagnostics completed {stage}")
                if status.get("status") == "failed":
                    failed.append(dict(label=row["label"], status=status))
            if all_done:
                journal(
                    root,
                    "intent",
                    "Write targeted relay statistics plots and recommendation packet",
                )
                make_report(prepared.parent.parent, cfg, output)
                journal(
                    root,
                    "outcome",
                    "Targeted relay packet complete; no PMF/convergence claim",
                )
                return
            if failed:
                raise RuntimeError(
                    f"Native target failed; partial accepted stages retained: {failed}"
                )
            if not workers_alive(prepared.parent.resolve()):
                raise RuntimeError(
                    "Incomplete native stages and no live target workers"
                )
            time.sleep(poll)
        raise TimeoutError(
            "24-hour bounded watcher deadline; incomplete stages retained"
        )
    except Exception as error:
        journal(root, "intent", "Write explicit targeted relay failure/partial packet")
        (output / "targeted_relay_partial.json").write_text(
            json.dumps(
                dict(
                    status="partial_or_failed",
                    error=f"{type(error).__name__}: {error}",
                    accepted_derived_stages=sorted(map(str, completed)),
                ),
                indent=2,
            )
            + "\n"
        )
        journal(
            root,
            "outcome",
            f"Targeted relay analysis partial/failure: {type(error).__name__}: {error}",
        )
        raise


def main() -> None:
    """Dispatch bounded diagnostics using explicitly supplied frozen inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("prepared", "output", "journal-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--max-hours", type=float, default=24)
    parser.add_argument("--poll-seconds", type=float, default=120)
    a = parser.parse_args()
    run(a.prepared, a.output, a.journal_root, a.max_hours, a.poll_seconds)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
