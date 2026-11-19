#!/usr/bin/env python3
"""Monitor coherent NEB bands, export VMD progress, and queue bounded candidate trials."""

from __future__ import annotations

import argparse
import fcntl
import json
import logging
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.neb_xyz_snapshot import snapshot_complete_band_from_xyz  # noqa: E402
from scripts.acylenzyme_candidate_tools import (  # noqa: E402
    export_band,
    launch_candidate,
    prepare_candidate,
    release_candidate,
)
from scripts.record_agent_action import append_action  # noqa: E402
from scripts.run_acylenzyme_neb import write_json  # noqa: E402

LOGGER = logging.getLogger(__name__)


def find_sander_pid(directory: Path, proc_root: Path = Path("/proc")) -> int | None:
    """Find the unique live serial sander process with exactly this working directory.

    Args:
        directory: Expected native NEB working directory.
        proc_root: Process filesystem root; injectable for tests.

    Returns:
        Matching PID, or None before the native stage starts.

    Raises:
        RuntimeError: If multiple matching writers exist.
    """
    found = []
    for process in proc_root.iterdir():
        if not process.name.isdigit():
            continue
        try:
            executable = (process / "exe").resolve(strict=True)
            cwd = (process / "cwd").resolve(strict=True)
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue  # Processes may exit while enumerating /proc.
        if executable.name in {"sander", "sander.LES"} and cwd == directory.resolve():
            found.append(int(process.name))
    if len(found) > 1:
        raise RuntimeError("Multiple sander writers in native NEB directory")
    return found[0] if found else None


def observe_candidates(
    state: dict, report_index: int, cycle: int, dips: tuple[int, ...]
) -> list[int]:
    """Select at most three distinct images persistent across two observed reports.

    Args:
        state: Mutable persistent monitor state.
        report_index: Unique complete native report index.
        cycle: Proven following native convergence cycle.
        dips: Strict interior local energy dips in this report.

    Returns:
        Eligible image numbers, excluding already reserved/launched images.
    """
    if report_index == state.get("last_report"):
        return []
    previous = set(state.get("previous_dips", []))
    launched = set(state.get("reserved_images", []))
    selected = sorted(previous.intersection(dips) - launched) if cycle >= 5 else []
    selected = selected[: max(0, 3 - len(launched))]
    state.update(last_report=report_index, previous_dips=list(dips))
    return selected


def audit(
    journal: Path, phase: str, action: str, files: list[Path] | None = None
) -> None:
    """Append a required durable action record before/after monitor mutations.

    Args:
        journal: Active append-only action journal.
        phase: Intent, outcome or note.
        action: Nonsecret action/result description.
        files: Explicit source/input/artifact paths to hash.
    """
    append_action(
        journal,
        "acylenzyme_candidate_monitor",
        phase,
        action,
        "tool",
        "Authorized bounded NEB monitoring and provisional trial releases",
        "Preserve inputs/logs; stop only the specifically identified trial if requested",
        files=[str(path) for path in files] if files else None,
    )


def monitor_once(config_path: Path, journal: Path) -> str:
    """Inspect one coherent band and perform due exports/trial reservations.

    Args:
        config_path: Campaign configuration.
        journal: Active journal; failed journal writes stop dependent mutations.

    Returns:
        Explicit monitor outcome for operational logging.
    """
    campaign = config_path.parent
    native = campaign / "03_neb"
    pid = find_sander_pid(native)
    if pid is None:
        status_path = campaign / "status.json"
        if status_path.exists():
            phase = json.loads(status_path.read_text())["phase"]
            if phase in {"stopped", "finished"}:
                return f"campaign terminal: {phase}"
        return "waiting for live native NEB writer"
    state_path = campaign / "candidate_monitor_state.json"
    state = (
        json.loads(state_path.read_text())
        if state_path.exists()
        else {
            "reserved_images": [],
            "previous_dips": [],
            "last_export_cycle": None,
        }
    )
    cfg = json.loads(config_path.read_text())
    audit(journal, "intent", "Read coherent native band snapshot", [config_path])
    with tempfile.TemporaryDirectory(prefix="acylenzyme_neb_band_") as directory:
        snapshot = Path(directory) / "band"
        try:
            report = snapshot_complete_band_from_xyz(
                native / "minimized.mdout",
                native / "nebpath.xyz",
                snapshot,
                cfg,
                16,
                writer_pid=pid,
            )
        except (FileNotFoundError, ValueError, RuntimeError) as error:
            audit(journal, "outcome", f"Snapshot deferred: {error}")
            return f"snapshot deferred: {error}"
        audit(
            journal,
            "outcome",
            f"Coherent report {report.report_index}, cycle {report.cycle}",
            [snapshot / "band_snapshot.json"],
        )
        if report.cycle is None:
            raise RuntimeError("Snapshot helper returned an unproven cycle")
        if report.cycle < state.get("last_cycle", 0):
            raise RuntimeError(
                "Native cycle regressed; campaign state must be reviewed"
            )
        selected = observe_candidates(
            state, report.report_index, report.cycle, report.provisional_minimum_images
        )
        state["last_cycle"] = report.cycle
        if (
            state["last_export_cycle"] is None
            or report.cycle - state["last_export_cycle"] >= 10
        ):
            audit(journal, "intent", f"Export coherent VMD band cycle {report.cycle}")
            generation = export_band(
                cfg, snapshot, campaign / "vmd_progress", report.cycle
            )
            state["last_export_cycle"] = report.cycle
            audit(
                journal,
                "outcome",
                f"Published {generation}",
                [generation / "MANIFEST.sha256"],
            )
        for image in selected:
            audit(
                journal, "intent", f"Reserve and launch provisional image {image} trial"
            )
            candidate = prepare_candidate(snapshot, campaign, image)
            state["reserved_images"].append(image)
            write_json(
                state_path, state
            )  # Reserve durably before any subprocess launch.
            try:
                launch_candidate(config_path, candidate, journal)
            except Exception as error:
                audit(
                    journal,
                    "outcome",
                    f"Candidate {image} launch failed: {error}",
                    [candidate / "candidate.json"],
                )
                raise
            audit(
                journal,
                "outcome",
                f"Launched provisional image {image} independent trial",
                [candidate / "candidate.json", candidate / "job/command.json"],
            )
        audit(journal, "intent", "Persist monitor observation state")
        write_json(state_path, state)
        audit(
            journal,
            "outcome",
            f"Observed cycle {report.cycle}; reserved {state['reserved_images']}",
            [state_path],
        )
    return f"observed cycle {report.cycle}; provisional dips {report.provisional_minimum_images}"


def main(argv: list[str] | None = None) -> int:
    """Run periodic campaign monitoring or one independently queued candidate release.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Zero on normal completion; unexpected errors propagate.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--journal", type=Path, default=ROOT / "output_agent_actions/20261021.jsonl"
    )
    parser.add_argument("--release", type=Path)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--interval-seconds", type=float, default=60)
    args = parser.parse_args(argv)
    if args.interval_seconds < 1:
        parser.error("--interval-seconds must be at least 1")
    args.config = args.config.resolve()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.release is not None:
        audit(
            args.journal,
            "intent",
            f"Unbiased candidate release {args.release}",
            [args.config, args.release / "candidate.ncrst"],
        )
        try:
            release_candidate(args.config, args.release)
        except Exception as error:
            audit(args.journal, "outcome", f"Candidate release failed: {error}")
            raise
        audit(
            args.journal,
            "outcome",
            f"Candidate trial finished {args.release}",
            [args.release / "trial_acceptance.json"],
        )
        return 0
    audit(args.journal, "intent", "Acquire exclusive candidate monitor lock")
    with (args.config.parent / "candidate_monitor.lock").open("a") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                "Another candidate monitor holds the campaign lock"
            ) from error
        audit(args.journal, "outcome", "Exclusive monitor lock acquired")
        while True:
            outcome = monitor_once(args.config, args.journal)
            LOGGER.info("%s", outcome)
            if args.once or outcome.startswith("campaign terminal:"):
                return 0
            time.sleep(args.interval_seconds)


if __name__ == "__main__":
    raise SystemExit(main())
