"""Immutable, statically partitioned umbrella controllers with explicit adoption."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import time
from typing import Any

from .jobs import launch_job

ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    """Return SHA-256 for a retained file."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            result.update(block)
    return result.hexdigest()


def read(path: Path) -> Any:
    """Read complete JSON, failing explicitly on partial or missing receipts."""
    return json.loads(path.read_text())


def write(path: Path, value: Any) -> None:
    """Atomically replace a controller status record."""
    pending = path.with_suffix(".pending")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def partition(schedule: list[str], running: list[str]) -> dict[str, list[list[dict]]]:
    """Assign eight adopted plus eight local and 47 remote windows exactly once."""
    if (
        len(schedule) != 63
        or len(set(schedule)) != 63
        or len(running) != 8
        or len(set(running)) != 8
        or not set(running) <= set(schedule)
    ):
        raise ValueError("Require 63 unique windows and eight unique running windows")
    remaining = [label for label in schedule if label not in running]
    local = [
        [{"label": a, "mode": "adopt"}, {"label": b, "mode": "launch"}]
        for a, b in zip(running, remaining[:8], strict=True)
    ]
    remote = [
        [{"label": label, "mode": "launch"} for label in remaining[8:][i::24]]
        for i in range(24)
    ]
    return {"local": local, "remote": remote}


def command(cfg: dict, prepared: Path, label: str) -> list[str]:
    """Return the original immutable window driver's exact argv."""
    return [
        cfg["python"],
        cfg["runtime_script"],
        "window",
        "--config",
        str(prepared),
        "--window",
        label,
    ]


def verify_job(cfg: dict, prepared: Path, label: str) -> dict | None:
    """Require existing adopted job identity to match the intended native window."""
    job = prepared.parent / label / "job"
    expected = command(cfg, prepared, label)
    metadata = read(job / "command.json")
    raw = (job / "status.json").read_bytes()
    status = json.loads(raw) if raw.endswith(b"\n") else None
    hashes = {str(prepared): digest(prepared), **cfg["prepared_input_sha256"]}
    for record in (metadata,) if status is None else (metadata, status):
        if record["command"] != expected or record["input_sha256"] != hashes:
            raise ValueError(f"Window job identity mismatch: {label}")
    return status


def complete(cfg: dict, prepared: Path, label: str) -> bool:
    """Distinguish native process success from a validated full-window receipt."""
    status = verify_job(cfg, prepared, label)
    if status is None or status["ended_utc"] is None:
        return False
    if status["exit_status"] != 0:
        raise RuntimeError(f"Window process failed: {label}")
    receipt = read(prepared.parent / label / "window_status.json")
    row = next(row for row in cfg["windows"] if row["label"] == label)
    if (
        receipt["status"] != "completed"
        or receipt["production_ps"] != 30
        or receipt["completed"] != [s["label"] for s in row["stages"]]
    ):
        raise ValueError(f"Incomplete scientific window receipt: {label}")
    return True


def validate_plan(plan: dict) -> None:
    """Reject duplicate, missing, or changed static ownership."""
    assignments = plan["assignments"]
    local, remote = plan["local_host"], plan["remote_host"]
    if (
        local == remote
        or set(assignments) != {local, remote}
        or len(assignments[local]) != 8
        or len(assignments[remote]) != 24
    ):
        raise ValueError("Require exactly eight local and 24 remote controllers")
    tasks = [
        task for workers in assignments.values() for tasks in workers for task in tasks
    ]
    labels = [task["label"] for task in tasks]
    if (
        len(labels) != 63
        or len(set(labels)) != 63
        or set(labels) != set(plan["windows"])
    ):
        raise ValueError("Static ownership must cover all 63 windows exactly once")
    for host, workers in assignments.items():
        for tasks in workers:
            modes = [task["mode"] for task in tasks]
            if (host == local and modes != ["adopt", "launch"]) or (
                host == remote and modes not in (["launch"], ["launch", "launch"])
            ):
                raise ValueError("Invalid static worker modes or capacity")


def prepare(master: Path, hosts: dict, local: str, remote: str) -> Path:
    """Freeze an explicit host plan and controller source without changing the campaign."""
    cfg = read(master)
    if digest(master) != (master.parent / "prepared.sha256").read_text().strip():
        raise ValueError("Master prepared digest mismatch")
    if set(hosts) != {local, remote} or local == remote:
        raise ValueError("Explicit distinct local/remote host identities required")
    if hosts[local] != {"prepared": str(master), "sha256": digest(master)}:
        raise ValueError("Local host identity must match master prepared config")
    running = [
        read(master.parent / "workers" / f"worker_{i:02d}" / "status.json")["current"]
        for i in range(1, 9)
    ]
    groups = partition(cfg["priority_schedule"], running)
    for label in running:
        if verify_job(cfg, master, label) is None:
            raise ValueError(f"Adopted window status write in progress: {label}")
    plan = {
        "hosts": hosts,
        "local_host": local,
        "remote_host": remote,
        "windows": cfg["priority_schedule"],
        "assignments": {local: groups["local"], remote: groups["remote"]},
    }
    validate_plan(plan)
    target = master.parent / "dispatch"
    target.mkdir(exist_ok=False)
    snapshot = target / "source_snapshot"
    names = [
        "scripts/dispatch_neb_umbrella.py",
        "src/nyl12_covalent/__init__.py",
        "src/nyl12_covalent/jobs.py",
        "src/nyl12_covalent/neb_umbrella_dispatch.py",
    ]
    for name in names:
        destination = snapshot / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, destination)
    plan["source_sha256"] = {str(p): digest(p) for p in snapshot.rglob("*.py")}
    plan["script"] = str(snapshot / names[0])
    path = target / "plan.json"
    write(path, plan)
    (target / "plan.sha256").write_text(digest(path) + "\n")
    return path


def frozen(path: Path, host: str) -> tuple[dict, dict, Path]:
    """Check host, plan, executable, original runtime, and all scientific input hashes."""
    if host != socket.gethostname().split(".")[0]:
        raise ValueError("Execution host does not match explicit plan host")
    if digest(path) != (path.parent / "plan.sha256").read_text().strip():
        raise ValueError("Dispatch plan changed")
    plan = read(path)
    validate_plan(plan)
    identity = plan["hosts"][host]
    prepared = Path(identity["prepared"])
    if digest(prepared) != identity["sha256"]:
        raise ValueError("Host prepared identity changed")
    if digest(prepared) != (prepared.parent / "prepared.sha256").read_text().strip():
        raise ValueError("Host prepared receipt changed")
    cfg = read(prepared)
    for name, expected in {
        **plan["source_sha256"],
        **cfg["frozen_source_sha256"],
        **cfg["execution_sha256"],
        **cfg["prepared_input_sha256"],
    }.items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Frozen source/input changed: {name}")
    if set(plan["windows"]) != {row["label"] for row in cfg["windows"]}:
        raise ValueError("Host window set differs from dispatch plan")
    return plan, cfg, prepared


def environment(cfg: dict) -> dict[str, str]:
    """Match the existing single-thread CPU window environment."""
    amber = Path(cfg["sander"]).parent.parent
    return {
        "AMBERHOME": str(amber),
        "LD_LIBRARY_PATH": str(amber / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "CUDA_VISIBLE_DEVICES": "",
    }


def worker(path: Path, host: str, index: int, poll: int = 30) -> None:
    """Adopt or launch each assigned window once, stopping on any failure."""
    plan, cfg, prepared = frozen(path, host)
    tasks = plan["assignments"][host]
    if not 0 <= index < len(tasks) or poll <= 0:
        raise ValueError("Invalid worker index or polling interval")
    folder = path.parent / host / f"worker{index + 1:02d}"
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "claim").open("x") as stream:
        stream.write(str(os.getpid()) + "\n")
    completed = []
    try:
        for task in tasks[index]:
            label = task["label"]
            frozen(path, host)
            write(
                folder / "status.json",
                {
                    "status": "running",
                    "current": label,
                    "mode": task["mode"],
                    "completed": completed,
                },
            )
            if task["mode"] == "launch":
                if (prepared.parent / label / "window_status.json").exists():
                    raise FileExistsError(f"Window already attempted: {label}")
                launch_job(
                    f"{cfg['session_prefix']}_{label}",
                    prepared.parent / label / "job",
                    prepared.parent,
                    command(cfg, prepared, label),
                    env_overrides=environment(cfg),
                    input_paths=[prepared, *map(Path, cfg["prepared_input_sha256"])],
                )
            while not complete(cfg, prepared, label):
                time.sleep(poll)
            completed.append(label)
        write(
            folder / "status.json",
            {
                "status": "completed",
                "completed": completed,
                "sampling_convergence": "Not asserted",
            },
        )
    except BaseException as error:
        write(
            folder / "status.json",
            {
                "status": "failed",
                "completed": completed,
                "error": f"{type(error).__name__}: {error}",
            },
        )
        raise


def launch(path: Path, host: str) -> None:
    """Launch only this host's statically assigned independent worker controllers."""
    plan, cfg, prepared = frozen(path, host)
    for index in range(len(plan["assignments"][host])):
        launch_job(
            f"{cfg['session_prefix']}_dispatch_{host}_{index + 1:02d}",
            path.parent / host / f"worker{index + 1:02d}" / "job",
            prepared.parent,
            [
                cfg["python"],
                plan["script"],
                "worker",
                "--plan",
                str(path),
                "--host",
                host,
                "--worker",
                str(index),
            ],
            env_overrides=environment(cfg),
            input_paths=[path, *map(Path, plan["source_sha256"])],
        )
