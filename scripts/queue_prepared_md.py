#!/usr/bin/env python3
"""Prepare/launch two sequential GPU queues of twelve matched100ns MD jobs."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.prepared_md_queue_plan import (  # noqa:E402
    environment,
    event,
    frozen_queue,
    plan_queue,
    runtime,
    write_json,
)


def worker(path: Path, repo: Path, gpu: str, poll: int = 30) -> None:
    """Run this GPU's tasks sequentially; each child has the existing tmux job receipt."""
    runner = runtime(repo)
    output, record = frozen_queue(path, repo)
    from nyl12_covalent.jobs import launch_job

    if gpu not in record["tasks"] or poll not in (30, 60):
        raise ValueError("Explicit planned GPU and30/60spoll required")
    folder = output / f"worker_gpu{gpu}"
    folder.mkdir(exist_ok=True)
    with (folder / "worker.claim").open("x") as stream:
        stream.write(str(os.getpid()) + "\n")
    completed = 0
    current = None
    try:
        runner.require_gpu(gpu)
        if os.environ.get("CUDA_VISIBLE_DEVICES") != gpu:
            raise ValueError(
                "Worker CUDA_VISIBLE_DEVICES differs from its explicit GPU"
            )
        for current in record["tasks"][gpu]:
            config = Path(current["config"])
            runner.frozen(config, repo)
            replica = Path(current["output"]) / current["replica"]
            job = replica / "job"
            write_json(
                folder / "status.json",
                {"status": "running", "completed": completed, "current": current},
            )
            event(folder, "replica_launch_intent", task=current)
            launch_job(
                current["session"],
                job,
                repo,
                [
                    sys.executable,
                    str(Path(runner.__file__).resolve()),
                    "execute",
                    "--repo",
                    str(repo),
                    "--config",
                    str(config),
                    "--replica",
                    current["replica"],
                ],
                env_overrides=environment(record, gpu),
                input_paths=[
                    output / "queue_plan.json",
                    *[Path(p) for p in record["source_sha256"]],
                ],
            )
            while True:
                raw = (job / "status.json").read_bytes()
                status = json.loads(raw) if raw.endswith(b"\n") else None
                if status is not None and status.get("ended_utc") is not None:
                    break
                time.sleep(poll)
            if status["exit_status"] != 0:
                raise RuntimeError(
                    f"Replica failed: {current['dataset']}/{current['replica']} exit{status['exit_status']}"
                )
            acceptance = json.loads(
                (replica / "100ns_process_complete.json").read_text()
            )
            if acceptance.get("completed") is not True:
                raise ValueError("Replica lacks successful100nsprocess receipt")
            completed += 1
            event(folder, "replica_completed", task=current, completed=completed)
        write_json(
            folder / "status.json",
            {
                "status": "completed",
                "completed": completed,
                "scientific_acceptance": "Not inferred from process completion",
            },
        )
    except BaseException as error:
        write_json(
            folder / "status.json",
            {
                "status": "failed",
                "completed": completed,
                "current": current,
                "error": f"{type(error).__name__}: {error}",
            },
        )
        event(
            folder,
            "worker_failed",
            task=current,
            error=f"{type(error).__name__}: {error}",
        )
        raise


def launch(path: Path, repo: Path) -> None:
    """Launch exactly two logged GPU workers; never call all-replica parallel launch."""
    runner = runtime(repo)
    output, record = frozen_queue(path, repo)
    from nyl12_covalent.jobs import launch_job

    for gpu in record["tasks"]:
        runner.require_gpu(gpu)
    for gpu in record["tasks"]:
        folder = output / f"worker_gpu{gpu}"
        folder.mkdir(exist_ok=True)
        launch_job(
            f"{record['config']['session_prefix']}_gpu{gpu}",
            folder / "job",
            repo,
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "worker",
                "--config",
                str(path),
                "--repo",
                str(repo),
                "--gpu",
                gpu,
            ],
            env_overrides=environment(record, gpu),
            input_paths=[
                output / "queue_plan.json",
                *[Path(p) for p in record["source_sha256"]],
            ],
        )
    write_json(
        output / "status.json",
        {"status": "two_workers_launched", "gpus": list(record["tasks"])},
    )


def main() -> None:
    """Expose explicit plan/launch/worker actions; preparation never launches jobs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "launch", "worker"))
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--existing-plans", action="store_true")
    parser.add_argument("--gpu")
    args = parser.parse_args()
    path = args.config.resolve(strict=True)
    repo = args.repo.resolve(strict=True)
    if args.action == "plan":
        plan_queue(path, repo, args.existing_plans)
    elif args.action == "launch":
        launch(path, repo)
    else:
        if args.gpu is None:
            parser.error("worker requires --gpu")
        worker(path, repo, args.gpu)


if __name__ == "__main__":
    main()
