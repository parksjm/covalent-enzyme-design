#!/usr/bin/env python3
"""Repair derived relay artifacts and resume a cleanly held native queue."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import importlib.util
import json
import subprocess
from pathlib import Path
import sys
import time


def journal(phase: str, action: str) -> None:
    """Record automated actions in the established append-only remote journal."""
    root = Path("/data/jpc/covalent-enzyme-design")
    subprocess.run(
        [
            "python3",
            str(root / "scripts/record_agent_action.py"),
            "--log",
            str(root / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_derived_guard",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Authorized derived-only repair and analysis-gated native continuation",
            "--rollback",
            "Preserve native artifacts and failure archives; no source hotedits",
        ],
        check=True,
    )


def load_queue(plan: Path) -> object:
    """Load the unchanged validated native controller for continuation only."""
    source = plan.parent / "source_snapshot/scripts/run_relay_umbrella.py"
    spec = importlib.util.spec_from_file_location("native_relay_queue", source)
    if spec is None or spec.loader is None:
        raise ValueError("Native queue source unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def continue_worker(plan: Path, output: Path, index: int) -> None:
    """Use original claims/inputs with fresh receipts and wait on derived holds."""
    queue = load_queue(plan)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "claim").open("x") as stream:
        stream.write(str(__import__("os").getpid()))
    while True:
        task = queue.claim_task(plan, index)
        if task is None:
            (output / "completed.json").write_text(
                json.dumps(dict(completed_utc=datetime.now(timezone.utc).isoformat()))
                + "\n"
            )
            return
        if isinstance(task, str):
            time.sleep(60)
            continue
        label, limit = task
        queue.event(output, "continuation_segment_intent", window=label, limit=limit)
        try:
            queue.run_segment(plan, label, limit)
        except Exception as error:
            queue.finish_claim(plan, label, limit, f"{type(error).__name__}: {error}")
            raise
        queue.finish_claim(plan, label, limit, None)
        queue.event(output, "continuation_segment_accepted", window=label, limit=limit)


def inspect_held(plan: Path) -> bool:
    """Require all original workers closed cleanly and all claims finalized."""
    cfg = json.loads(plan.read_text())
    statuses = list((plan.parent / "workers").glob("*/job/status.json"))
    if len(statuses) != cfg["workers"]:
        if (
            "initial_workers" in cfg
            and cfg["initial_workers"] <= len(statuses) < cfg["workers"]
        ):
            return False
        raise ValueError("Original worker inventory incomplete")
    jobs = [json.loads(p.read_text()) for p in statuses]
    if any(row["ended_utc"] is None for row in jobs):
        return False
    if any(row["exit_status"] != 0 for row in jobs):
        raise RuntimeError("Native worker failure requires review")
    claims = [
        json.loads(p.read_text()) for p in (plan.parent / "claims").glob("*.json")
    ]
    if any(row["status"] != "completed" for row in claims):
        raise RuntimeError("Unresolved native queue claim requires review")
    states = [
        json.loads(p.read_text())
        for p in (plan.parent / "workers").glob("*/status.json")
    ]
    return any(row["status"] == "held_for_analysis" for row in states)


def guard(args: argparse.Namespace) -> None:
    """Repair accepted outputs, optionally without inspecting the native queue."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    import parmed
    from nyl12_covalent.jobs import launch_job
    from nyl12_covalent.relay_umbrella_runtime import ensure_diagnostics
    from nyl12_covalent.relay_umbrella_analysis_io import load_block

    cfg = json.loads((args.campaign / "prepared.json").read_text())
    plan = json.loads(args.plan.read_text())
    request = json.loads((args.campaign / "analysis_request.json").read_text())
    top = parmed.load_file(str(args.campaign / cfg["topology"]))
    args.output.mkdir(exist_ok=False)
    with (args.plan.parent / "derived_guard.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        deadline = time.monotonic() + args.duration_hours * 3600
        while time.monotonic() < deadline:
            results = []
            for row in cfg["windows"]:
                for stage in row["stages"]:
                    if stage["kind"] != "production":
                        continue
                    folder = args.campaign / Path(stage["input"]).parent
                    if not (folder / "acceptance.json").is_file():
                        continue
                    outcome = folder / "analysis_outcome.json"
                    if (
                        outcome.exists()
                        and json.loads(outcome.read_text())["status"] == "completed"
                    ):
                        try:
                            load_block(
                                folder / "analysis_acceptance.json",
                                temperature=cfg["temperature_K"],
                                hamiltonian=cfg["physical_hamiltonian_sha256"],
                                diagnostic_columns=request["diagnostic_columns"],
                            )
                            continue
                        except (KeyError, ValueError, FileNotFoundError) as error:
                            journal(
                                "intent",
                                "Preserve and invalidate legacy derived receipt "
                                + str(folder)
                                + " "
                                + str(error),
                            )
                            archive = folder / (
                                "analysis_outcome_before_revalidation_"
                                + datetime.now(timezone.utc).strftime(
                                    "%Y%m%dT%H%M%S%fZ"
                                )
                                + ".json"
                            )
                            archive.write_bytes(outcome.read_bytes())
                            outcome.write_text(
                                json.dumps(
                                    dict(
                                        status="failed",
                                        error="Strict native-parity loader revalidation: "
                                        + str(error),
                                    )
                                )
                                + "\n"
                            )
                            journal(
                                "outcome",
                                "Legacy invalid derived outcome preserved; native unchanged",
                            )
                    # Original native worker must finish its own postprocess attempt first.
                    receipt = (
                        args.campaign
                        / row["label"]
                        / f"segment_receipt_limit_{sum(s['kind'] == 'production' for s in row['stages'][: row['stages'].index(stage) + 1]):02d}.json"
                    )
                    if not receipt.is_file():
                        continue
                    journal(
                        "intent",
                        "Regenerate derived-only relay artifacts " + str(folder),
                    )
                    result = ensure_diagnostics(cfg, row, folder, top)
                    journal(
                        "outcome",
                        "Derived-only repair completed="
                        + str(result)
                        + " "
                        + str(folder),
                    )
                    results.append(
                        dict(
                            window=row["label"],
                            stage=stage["label"],
                            analysis_completed=result,
                        )
                    )
                    if result:
                        load_block(
                            folder / "analysis_acceptance.json",
                            temperature=cfg["temperature_K"],
                            hamiltonian=cfg["physical_hamiltonian_sha256"],
                            diagnostic_columns=request["diagnostic_columns"],
                        )
                    if not result:
                        raise RuntimeError(f"Derived repair failed: {folder}")
            (args.output / "latest.json").write_text(
                json.dumps(
                    dict(
                        observed_utc=datetime.now(timezone.utc).isoformat(),
                        repairs=results,
                    ),
                    indent=2,
                )
                + "\n"
            )
            if args.analysis_only:
                time.sleep(60)
                continue
            if (
                inspect_held(args.plan)
                and not (args.output / "continuations_launched.json").exists()
            ):
                # Recheck every first-block analysis before adding any native task.
                all_ready = True
                for row in cfg["windows"]:
                    first = next(s for s in row["stages"] if s["kind"] == "production")
                    path = (
                        args.campaign
                        / Path(first["input"]).parent
                        / "analysis_outcome.json"
                    )
                    all_ready &= (
                        path.is_file()
                        and json.loads(path.read_text())["status"] == "completed"
                    )
                if all_ready:
                    commands = []
                    journal(
                        "intent",
                        "Launch cleanly-held native queue continuation with "
                        + str(plan["workers"])
                        + " workers",
                    )
                    for index in range(plan["workers"]):
                        destination = args.output / f"continuation_{index:02d}"
                        argv = [
                            plan["python"],
                            str(Path(__file__).resolve()),
                            "worker",
                            "--plan",
                            str(args.plan),
                            "--output",
                            str(destination / "receipts"),
                            "--index",
                            str(index),
                        ]
                        launch_job(
                            f"relay_us_continuation_{index:02d}_20261125",
                            destination / "job",
                            args.campaign,
                            argv,
                            env_overrides=plan["environment"],
                            input_paths=[args.plan, Path(__file__)],
                        )
                        commands.append(argv)
                    journal("outcome", "Native queue continuation workers launched")
                    (args.output / "continuations_launched.json").write_text(
                        json.dumps(dict(argv=commands), indent=2) + "\n"
                    )
            time.sleep(60)


def main() -> None:
    """Run an explicitly bounded guardian or one continuation worker."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("guard", "worker"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--campaign", type=Path)
    parser.add_argument("--duration-hours", type=float)
    parser.add_argument("--index", type=int)
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Repair derived receipts without inspecting or launching native workers",
    )
    args = parser.parse_args()
    if args.action == "worker":
        if args.analysis_only:
            parser.error("analysis-only is available only for guard")
        if args.index is None:
            parser.error("worker requires index")
        continue_worker(args.plan, args.output, args.index)
    else:
        if (
            args.campaign is None
            or args.duration_hours is None
            or args.duration_hours <= 0
        ):
            parser.error("guard requires campaign and positive bounded duration")
        guard(args)


if __name__ == "__main__":
    main()
