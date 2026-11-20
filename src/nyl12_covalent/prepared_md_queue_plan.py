#!/usr/bin/env python3
"""Prepare/launch two sequential GPU queues of twelve matched100ns MD jobs."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
from typing import Any


def runtime(repo: Path) -> Any:
    """Load the supplied repository's portable runner without copying MD logic."""
    sys.path.insert(0, str(repo / "student_handoff/nyl12_ash274_md/tools"))
    module = importlib.import_module("run_prepared_md")
    module.context(repo)
    return module


def digest(path: Path) -> str:
    """Hash a small queue/plan source; large scientific inputs use runner hashes."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def write_json(path: Path, value: Any) -> None:
    """Atomically publish a mutable queue status or new frozen plan."""
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def event(folder: Path, action: str, **fields: Any) -> None:
    """Append one durable worker intent/outcome receipt."""
    with (folder / "events.jsonl").open("a") as stream:
        stream.write(
            json.dumps(
                {
                    "utc": datetime.now(timezone.utc).isoformat(),
                    "event": action,
                    **fields,
                }
            )
            + "\n"
        )
        stream.flush()
        os.fsync(stream.fileno())


def plan_queue(path: Path, repo: Path, existing_plans: bool = False) -> Path:
    """Freeze three conditions × four replicas into balanced per-GPU queues.

    Args:
        path: JSON with configs[3],gpus[2],output and session_prefix.
        repo: Explicit code checkout containing the portable prepared-MD tools.
        existing_plans: Require existing validated MD plans instead of creating new.

    Returns:
        New queue directory; no tmux or MD process is launched.

    Raises:
        ValueError: If conditions, seeds, GPU allocation or engines disagree.
    """
    runner = runtime(repo)
    cfg = json.loads(path.read_text())
    gpus = list(map(str, cfg["gpus"]))
    if (
        len(cfg["configs"]) != 3
        or len(gpus) != 2
        or len(set(gpus)) != 2
        or any(not g.isdecimal() for g in gpus)
    ):
        raise ValueError(
            "Require three condition configs and two distinct numeric GPU IDs"
        )
    entrypoint = Path(__file__).resolve().parents[2] / "scripts/queue_prepared_md.py"
    groups = {g: [] for g in gpus}
    hashes = {
        str(path): digest(path),
        str(Path(__file__).resolve()): digest(Path(__file__)),
    }
    hashes[str(entrypoint)] = digest(entrypoint)
    conditions = set()
    used_seeds = set()
    environments = set()
    for item in cfg["configs"]:
        config = (path.parent / item).resolve(strict=True)
        if not existing_plans:
            runner.plan(config, repo)
        md, output, record = runner.frozen(config, repo)
        runner.replicas(md)
        if md["dataset"] in conditions or len(md["replicas"]) != 4:
            raise ValueError(
                "Require three distinct conditions, each with four replicas"
            )
        conditions.add(md["dataset"])
        home, engine = runner.engines(md, config.parent)
        environments.add((str(home), str(engine)))
        hashes.update(record["source_sha256"])
        hashes[str(output / "plan.json")] = digest(output / "plan.json")
        for gpu in gpus:
            rows = [r for r in md["replicas"] if str(r["gpu"]) == gpu]
            if len(rows) != 2:
                raise ValueError(
                    "Each condition must allocate exactly two replicas to each GPU"
                )
            tasks = []
            for row in rows:
                seeds = {row["seed"] + 10000 * i for i in range(6)}
                if seeds & used_seeds:
                    raise ValueError("Stage seeds overlap across conditions/replicas")
                used_seeds.update(seeds)
                tasks.append(
                    {
                        "config": str(config),
                        "output": str(output),
                        "replica": row["name"],
                        "dataset": md["dataset"],
                        "seed": row["seed"],
                        "gpu": gpu,
                        "session": f"{md['session_prefix']}_{row['name']}",
                    }
                )
            groups[gpu].append(tasks)
    if len(environments) != 1:
        raise ValueError(
            "Matched conditions require the same explicit Amber/CUDA installation"
        )
    tasks = {}
    for index, gpu in enumerate(gpus):
        order = groups[gpu][index:] + groups[gpu][:index]
        tasks[gpu] = [condition[replica] for replica in range(2) for condition in order]
    output = (path.parent / cfg["output"]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    home, engine = next(iter(environments))
    record = {
        "config": cfg,
        "queue_config": str(path),
        "repo": str(repo),
        "tasks": tasks,
        "source_sha256": hashes,
        "amber_home": home,
        "cuda_engine": engine,
        "policy": "One serial worker per GPU; fail that worker on any job/scientific error; other GPU independent",
        "production_ns_per_job": 100,
        "total_jobs": 12,
        "unique_stage_seeds": len(used_seeds),
    }
    write_json(output / "queue_plan.json", record)
    (output / "queue_plan.sha256").write_text(digest(output / "queue_plan.json") + "\n")
    write_json(
        output / "status.json", {"status": "prepared_not_launched", "total_jobs": 12}
    )
    return output


def frozen_queue(path: Path, repo: Path) -> tuple[Path, dict[str, Any]]:
    """Require unchanged queue config, accepted plans and all pinned source inputs."""
    cfg = json.loads(path.read_text())
    output = (path.parent / cfg["output"]).resolve(strict=True)
    plan = output / "queue_plan.json"
    if digest(plan) != (output / "queue_plan.sha256").read_text().strip():
        raise ValueError("Frozen queue plan changed")
    record = json.loads(plan.read_text())
    if record["config"] != cfg or record["repo"] != str(repo):
        raise ValueError("Queue config/repository differs from frozen plan")
    for source, expected in record["source_sha256"].items():
        if digest(Path(source)) != expected:
            raise ValueError(f"Frozen queue source/input changed: {source}")
    return output, record


def environment(record: dict[str, Any], gpu: str) -> dict[str, str]:
    """Return the same explicit engine environment for worker and replica jobs."""
    return {
        "AMBERHOME": record["amber_home"],
        "LD_LIBRARY_PATH": str(Path(record["amber_home"]) / "lib"),
        "CUDA_VISIBLE_DEVICES": gpu,
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
