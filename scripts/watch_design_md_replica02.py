"""Prepare and run a source-pinned serial queue for completed replica02 MD."""

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

LOGGER = logging.getLogger(__name__)


class JournalError(RuntimeError):
    """Stop dependent work if the mandatory action journal fails."""


FOLLOW = Path("design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup")
DATASETS = {"WT", "broader_mask07", "focused_mask03", "V111C"}


def digest(path: Path) -> str:
    """Return SHA-256 of a file without loading it into memory."""
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Atomically publish a JSON receipt."""
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(temp, path)


def prepare(
    project: Path,
    queue: Path,
    python: Path,
    cases: list[str],
    sources: list[Path],
    code_root: Path,
    legacy: Path,
) -> dict[str, Any]:
    """Create a host queue, preserving exact snapshots and source digests.

    Args:
        project: Absolute repository root on the target host.
        queue: New host queue directory.
        python: Explicit host Python executable.
        cases: Exactly two distinct assigned dataset names.
        sources: Explicit absolute analyzer and dependency source paths.
        code_root: Absolute repository-shaped source tree used by the analyzer.
        legacy: Absolute existing tier-classifier source file.

    Returns:
        The retained source-pinned host plan.
    """
    if not project.is_absolute() or not python.is_absolute() or not python.is_file():
        raise ValueError("Absolute project and existing Python are required")
    if len(cases) != 2 or len(set(cases)) != 2 or not set(cases) <= DATASETS:
        raise ValueError("Require two distinct supported datasets")
    if (
        not sources
        or len(set(sources)) != len(sources)
        or any(not p.is_absolute() or not p.is_file() for p in sources)
    ):
        raise ValueError("Require distinct absolute existing source paths")
    if (
        not code_root.is_absolute()
        or not legacy.is_absolute()
        or legacy.name != "analyze_nyl12_nma_four_site_geometry.py"
    ):
        raise ValueError("Explicit code root and legacy classifier required")
    required = {
        code_root / "scripts/analyze_design_md_geometry.py",
        code_root / "scripts/watch_design_md_replica02.py",
        code_root / "src/nyl12_covalent/__init__.py",
        code_root / "src/nyl12_covalent/design_md_geometry.py",
        code_root / "src/nyl12_covalent/design_md_completion.py",
        code_root / "src/nyl12_covalent/active_site_geometry.py",
        code_root / "src/nyl12_covalent/catalytic_alignment.py",
        code_root / "src/nyl12_covalent/catalytic_alignment_geometry.py",
        code_root / "src/nyl12_covalent/native_thr_policy.py",
        code_root
        / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation_source_snapshot_checked/scripts/screen_pa6_binding_orientation.py",
        code_root / "scripts/record_agent_action.py",
        code_root / "scripts/run_nyl12_job.py",
        code_root / "src/nyl12_covalent/jobs.py",
        legacy,
    }
    if not required <= set(sources):
        raise ValueError(
            f"Missing required source files: {sorted(map(str, required - set(sources)))}"
        )
    queue.mkdir(parents=True, exist_ok=False)
    snap = queue / "source_snapshot"
    snap.mkdir()
    records = []
    for index, source in enumerate(sources):
        target = snap / f"{index:03d}_{source.name}"
        expected = digest(source)
        shutil.copy2(source, target)
        if digest(source) != expected or digest(target) != expected:
            raise ValueError(f"Source changed during snapshot: {source}")
        records.append(
            {
                "path": str(source),
                "snapshot": str(target.relative_to(queue)),
                "sha256": expected,
            }
        )
    plan = {
        "project": str(project),
        "code_root": str(code_root),
        "legacy": str(legacy),
        "python": str(python),
        "python_sha256": digest(python),
        "cases": cases,
        "sources": records,
        "replica": "02",
        "native_thr_policy": "diagnostic_only",
    }
    write_json(queue / "plan.json", plan)
    return plan


def verify_sources(queue: Path, plan: dict[str, Any]) -> None:
    """Reject any changed live source or retained snapshot."""
    if digest(Path(plan["python"])) != plan["python_sha256"]:
        raise ValueError("Pinned Python executable changed")
    for record in plan["sources"]:
        for path in (Path(record["path"]), queue / record["snapshot"]):
            if digest(path) != record["sha256"]:
                raise ValueError(f"Pinned source changed: {path}")


def readiness(project: Path, dataset: str) -> str:
    """Distinguish pending, completed and failed native MD processes."""
    run = project / FOLLOW / "md" / dataset / "replica_02"
    marker = run / "100ns_process_complete.json"
    if marker.exists():
        if json.loads(marker.read_text()).get("completed") is not True:
            raise ValueError(f"Invalid completed MD receipt: {marker}")
        return "ready"
    status = run / "job/status.json"
    if status.exists():
        record = json.loads(status.read_text())
        if record.get("ended_utc") is not None and record.get("exit_status") != 0:
            return "failed"
    return "pending"


def _journal(
    project: Path, code_root: Path, phase: str, action: str, files: tuple[Path, ...]
) -> None:
    command = [
        sys.executable,
        str(code_root / "scripts/record_agent_action.py"),
        "--log",
        str(project / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "replica02_analysis_watcher",
        "--phase",
        phase,
        "--action",
        action,
        "--permission-mode",
        "escalated",
        "--reason",
        "PI-authorized completed-replica analysis",
        "--rollback",
        "Preserve MD inputs and case analysis evidence",
    ]
    for file in files:
        command.extend(("--file", str(file)))
    try:
        subprocess.run(command, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise JournalError("Mandatory action journal failed") from error


def run_case(queue: Path, plan: dict[str, Any], dataset: str) -> str:
    """Run exactly one completed case, retaining separate process and scientific records."""
    project = Path(plan["project"])
    code_root = Path(plan["code_root"])
    plan_identity = digest(queue / "plan.json")
    case = queue / dataset
    case.mkdir(exist_ok=False)
    output = project / FOLLOW / "analysis/replica_02" / dataset
    receipt = (
        project / FOLLOW / "md" / dataset / "replica_02/100ns_process_complete.json"
    )
    try:
        if output.exists():
            raise FileExistsError(output)
        verify_sources(queue, plan)
        identity = digest(receipt)
        command = [
            plan["python"],
            str(code_root / "scripts/analyze_design_md_geometry.py"),
            "--dataset",
            dataset,
            "--replica",
            "02",
            "--project-root",
            str(project),
            "--output",
            str(output),
        ]
        env = os.environ.copy()
        env.update(
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
            PYTHONPATH=f"{code_root}/src:{code_root}:{Path(plan['legacy']).parent}",
        )
        overrides = {
            key: env[key]
            for key in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "PYTHONPATH",
            )
        }
        write_json(
            case / "process.json",
            {
                "status": "running",
                "command": command,
                "environment_overrides": overrides,
                "completion_sha256": identity,
            },
        )
        _journal(
            project,
            code_root,
            "intent",
            f"Run completed-replica analysis: {command}",
            (queue / "plan.json", receipt),
        )
        with (case / "analysis.log").open("x") as log:
            process = subprocess.run(
                command,
                cwd=project,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        _journal(
            project,
            code_root,
            "outcome",
            f"Completed-replica analyzer exited {process.returncode}: {command}",
            (queue / "plan.json", receipt),
        )
        write_json(
            case / "process.json",
            {
                "status": "completed",
                "exit_status": process.returncode,
                "command": command,
                "environment_overrides": overrides,
                "completion_sha256": identity,
            },
        )
        verify_sources(queue, plan)
        if process.returncode != 0:
            raise RuntimeError(f"Analyzer exited {process.returncode}")
        if digest(receipt) != identity or digest(queue / "plan.json") != plan_identity:
            raise ValueError(
                "MD completion receipt or queue plan changed during analysis"
            )
        names = ("summary.json", "sites.csv", "sites.jsonl")
        if any(
            not (output / name).is_file() or not (output / name).stat().st_size
            for name in names
        ):
            raise ValueError("Analyzer outputs missing or empty")
        summary = json.loads((output / "summary.json").read_text())
        if (
            summary.get("dataset") != dataset
            or summary.get("replica") != "replica_02"
            or summary.get("native_thr_role") != "diagnostic_only"
            or summary.get("selection_policy")
            != "Tier A/B catalytic geometry only; native and joint metrics retained as historical diagnostics"
        ):
            raise ValueError("Analyzer output identity or native Thr policy mismatch")
        write_json(
            case / "scientific_acceptance.json",
            {
                "status": "accepted",
                "dataset": dataset,
                "replica": "replica_02",
                "native_thr_policy": "diagnostic_only",
                "completion_sha256": identity,
                "output_sha256": {name: digest(output / name) for name in names},
            },
        )
        return "accepted"
    except JournalError:
        raise
    except Exception as error:
        write_json(
            case / "scientific_acceptance.json",
            {"status": "failed", "dataset": dataset, "error": repr(error)},
        )
        LOGGER.exception("Analysis case %s failed", dataset)
        return "failed"


def watch(queue: Path, poll_seconds: float, timeout_hours: float) -> dict[str, str]:
    """Wait for assigned MD completions and run one CPU analysis at a time.

    Args:
        queue: Prepared host queue directory.
        poll_seconds: Positive polling period.
        timeout_hours: Positive total waiting limit.

    Returns:
        Terminal per-case status mapping.
    """
    if poll_seconds <= 0 or timeout_hours <= 0:
        raise ValueError("Polling interval and timeout must be positive")
    plan_path = queue / "plan.json"
    plan_hash = digest(plan_path)
    plan = json.loads(plan_path.read_text())
    verify_sources(queue, plan)
    state = {case: "pending" for case in plan["cases"]}
    deadline = time.monotonic() + timeout_hours * 3600
    while "pending" in state.values():
        for case in plan["cases"]:
            if state[case] != "pending":
                continue
            try:
                if digest(plan_path) != plan_hash:
                    raise ValueError("Queue plan changed")
                status = readiness(Path(plan["project"]), case)
                if status == "pending":
                    continue
                if status == "failed":
                    raise ValueError("Native MD job failed before completion")
                state[case] = run_case(queue, plan, case)
            except JournalError:
                raise
            except Exception as error:
                folder = queue / case
                folder.mkdir(exist_ok=False)
                write_json(
                    folder / "scientific_acceptance.json",
                    {"status": "failed", "dataset": case, "error": repr(error)},
                )
                state[case] = "failed"
        write_json(queue / "queue_status.json", {"cases": state})
        if "pending" not in state.values():
            break
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            state = {
                key: "timeout" if value == "pending" else value
                for key, value in state.items()
            }
            write_json(queue / "queue_status.json", {"cases": state})
            break
        time.sleep(min(poll_seconds, remaining))
    return state


def main() -> None:
    """Prepare or execute an explicitly configured replica02 host queue."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    setup = sub.add_parser("prepare")
    for option in ("project", "queue", "python"):
        setup.add_argument(f"--{option}", type=Path, required=True)
    setup.add_argument("--code-root", type=Path, required=True)
    setup.add_argument("--legacy", type=Path, required=True)
    setup.add_argument("--case", action="append", required=True)
    setup.add_argument("--source", action="append", type=Path, required=True)
    run = sub.add_parser("watch")
    run.add_argument("--queue", type=Path, required=True)
    run.add_argument("--poll-seconds", type=float, default=120)
    run.add_argument("--timeout-hours", type=float, default=48)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.action == "prepare":
        prepare(
            args.project,
            args.queue,
            args.python,
            args.case,
            args.source,
            args.code_root,
            args.legacy,
        )
    else:
        result = watch(args.queue, args.poll_seconds, args.timeout_hours)
        if any(value != "accepted" for value in result.values()):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
