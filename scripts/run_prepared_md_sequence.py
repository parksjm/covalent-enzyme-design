#!/usr/bin/env python3
"""Run an explicit sequence of existing prepared MD jobs on one reserved GPU."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.jobs import launch_job  # noqa: E402
from nyl12_covalent.prepared_md_queue_plan import (  # noqa: E402
    digest,
    runtime,
    write_json,
)

LOGGER = logging.getLogger(__name__)


def checked_tasks(path: Path, repo: Path, gpu: str) -> tuple[Any, list[dict[str, Any]]]:
    """Require frozen inputs, unique replicas/seeds and one explicit GPU.

    Args:
        path: Queue JSON containing tasks and a source_sha256 dictionary.
        repo: Pinned checkout or private source snapshot.
        gpu: NVIDIA device index reserved for this queue.

    Returns:
        Existing prepared-MD runtime and fully validated tasks.
    """
    document = json.loads(path.read_text())
    if document["gpu"] != gpu or os.environ.get("CUDA_VISIBLE_DEVICES") != gpu:
        raise ValueError("Queue, environment and requested GPU must agree")
    for source, expected in document["source_sha256"].items():
        if digest(Path(source)) != expected:
            raise ValueError(f"Queue source changed: {source}")
    runner = runtime(repo)
    tasks, identities, seeds = [], set(), set()
    for task in document["tasks"]:
        config_path = Path(task["config"]).resolve(strict=True)
        config, output, record = runner.frozen(config_path, repo)
        runner.replicas(config)
        row = next(r for r in config["replicas"] if r["name"] == task["replica"])
        identity = (str(output), row["name"])
        stage_seeds = {row["seed"] + 10000 * i for i in range(6)}
        if str(row["gpu"]) != gpu or identity in identities or seeds & stage_seeds:
            raise ValueError("Queue has a duplicate replica/seed or wrong GPU")
        identities.add(identity)
        seeds.update(stage_seeds)
        tasks.append(dict(config=config_path, row=row, output=output, record=record))
    if not tasks:
        raise ValueError("Queue must contain at least one explicit task")
    return runner, tasks


def require_idle_gpu(gpu: str) -> None:
    """Reject occupied GPUs instead of competing with another calculation."""
    result = subprocess.run(
        ["nvidia-smi", "-i", gpu, "--query-compute-apps=pid", "--format=csv,noheader"],
        check=True,
        capture_output=True,
        text=True,
    )
    if result.stdout.strip():
        raise RuntimeError(f"GPU {gpu} has an existing compute process")


def run_queue(path: Path, repo: Path, gpu: str, status_path: Path) -> None:
    """Launch each existing replica in tmux and advance only after acceptance.

    Args:
        path: Hashed queue configuration.
        repo: Code snapshot used by the original prepared MD entry point.
        gpu: Device index, also required in CUDA_VISIBLE_DEVICES.
        status_path: Mutable queue receipt, distinct from child receipts.
    """
    runner, tasks = checked_tasks(path, repo, gpu)
    completed = 0
    try:
        for task in tasks:
            require_idle_gpu(gpu)
            runner.require_gpu(gpu)
            config = task["config"]
            cfg, _, record = runner.frozen(config, repo)
            row = task["row"]
            folder = task["output"] / row["name"]
            job = folder / "job"
            write_json(
                status_path,
                dict(
                    status="running",
                    completed=completed,
                    dataset=cfg["dataset"],
                    replica=row["name"],
                ),
            )
            LOGGER.info("Starting %s/%s on GPU %s", cfg["dataset"], row["name"], gpu)
            launch_job(
                f"{cfg['session_prefix']}_{row['name']}",
                job,
                repo,
                [
                    sys.executable,
                    str(Path(runner.__file__)),
                    "execute",
                    "--repo",
                    str(repo),
                    "--config",
                    str(config),
                    "--replica",
                    row["name"],
                ],
                env_overrides={
                    key: os.environ[key]
                    for key in (
                        "AMBERHOME",
                        "LD_LIBRARY_PATH",
                        "CUDA_VISIBLE_DEVICES",
                        "OMP_NUM_THREADS",
                        "OPENBLAS_NUM_THREADS",
                        "MKL_NUM_THREADS",
                    )
                },
                input_paths=[path, *map(Path, record["source_sha256"])],
            )
            while True:
                state = json.loads((job / "status.json").read_text())
                if state["ended_utc"] is not None:
                    break
                time.sleep(60)
            if state["exit_status"] != 0:
                raise RuntimeError(f"MD child failed: {job}")
            receipt = json.loads((folder / "100ns_process_complete.json").read_text())
            if receipt["completed"] is not True:
                raise ValueError(f"Missing accepted 100 ns completion: {folder}")
            completed += 1
        write_json(
            status_path,
            dict(
                status="completed",
                completed=completed,
                scientific_acceptance="Not inferred from completion",
            ),
        )
    except BaseException as exc:
        write_json(
            status_path,
            dict(
                status="failed",
                completed=completed,
                error=f"{type(exc).__name__}: {exc}",
            ),
        )
        raise


def main() -> None:
    """Run only the explicitly configured, frozen queue."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--gpu", required=True)
    parser.add_argument("--status", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    run_queue(
        args.config.resolve(strict=True),
        args.repo.resolve(strict=True),
        args.gpu,
        args.status,
    )


if __name__ == "__main__":
    main()
