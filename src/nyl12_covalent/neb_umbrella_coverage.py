"""Schedule first-block coverage without changing native umbrella inputs."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import shlex
import socket
import subprocess
import time
from typing import Any

from .jobs import launch_job
from .neb_umbrella_campaign import digest, frozen, write_json
from .neb_umbrella_handoff import capture, signal_verified, suspend_driver

LOGGER = logging.getLogger(__name__)


def context(schedule: Path, host: str) -> tuple[dict, dict, Path]:
    """Validate scheduling sources and the unchanged host scientific preparation."""
    plan = json.loads(schedule.read_text())
    if digest(schedule) != schedule.with_suffix(".sha256").read_text().strip():
        raise ValueError("Coverage scheduling plan changed")
    if host != socket.gethostname().split(".")[0]:
        raise ValueError("Explicit host does not match execution host")
    for name, expected in plan["source_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Coverage source changed: {name}")
    prepared = Path(plan["hosts"][host]["prepared"])
    if digest(prepared) != plan["hosts"][host]["sha256"]:
        raise ValueError("Original scientific preparation changed")
    return plan, frozen(prepared), prepared


def event(folder: Path, action: str, **details: Any) -> None:
    """Append durable controller actions before and after process changes."""
    record = dict(utc=datetime.now(timezone.utc).isoformat(), action=action, **details)
    with (folder / "events.jsonl").open("a") as stream:
        stream.write(json.dumps(record, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def inventory(schedule: Path, host: str) -> Path:
    """Capture exact old controllers and their live native children without signaling."""
    plan, cfg, prepared = context(schedule, host)
    original = Path(plan["parent_plans"][host])
    parent = json.loads(original.read_text())
    native = {}
    for path in Path("/proc").iterdir():
        if not path.name.isdigit():
            continue
        try:
            argv = (path / "cmdline").read_bytes().split(b"\0")
            if not argv or argv[0].decode(errors="replace") != cfg["sander"]:
                continue
            item = capture(int(path.name))
        except (FileNotFoundError, ProcessLookupError, PermissionError) as error:
            LOGGER.debug("Unrelated/transient proc entry: %s", error)
            continue
        cwd = Path(item["cwd"])
        if cwd.parent.parent != prepared.parent:
            continue
        label = cwd.parent.name
        if label in native:
            raise ValueError(f"Multiple native processes for {label}")
        native[label] = item
    rows = []
    for index, tasks in enumerate(plan["assignments"][host]):
        label = tasks[0]["label"]
        current = native[label]
        driver = capture(current["ppid"])
        expected = [
            cfg["python"],
            cfg["runtime_script"],
            "window",
            "--config",
            str(prepared),
            "--window",
            label,
        ]
        if driver["argv"] != expected or driver["cwd"] != str(prepared.parent):
            raise ValueError(f"Native parent identity mismatch: {label}")
        claim = original.parent / host / f"worker{index + 1:02d}" / "claim"
        controller = capture(int(claim.read_text()))
        expected = [
            cfg["python"],
            parent["script"],
            "worker",
            "--plan",
            str(original),
            "--host",
            host,
            "--worker",
            str(index),
        ]
        if controller["argv"] != expected:
            raise ValueError(f"Old queue controller identity mismatch: {claim}")
        for task in tasks[1:]:
            folder = prepared.parent / task["label"]
            if (folder / "window_status.json").exists() or list(
                folder.glob("*/process.json")
            ):
                raise ValueError(f"Queued window already attempted: {folder}")
        rows.append(
            dict(
                current=label,
                queued=[t["label"] for t in tasks[1:]],
                native=current,
                driver=driver,
                controller=controller,
            )
        )
    if len(native) != len(rows):
        raise ValueError("Unexpected unowned native campaign process")
    folder = schedule.parent / host
    folder.mkdir(exist_ok=True)
    destination = folder / "inventory.json"
    if destination.exists():
        raise FileExistsError(destination)
    write_json(
        destination,
        dict(
            host=host,
            schedule_sha256=digest(schedule),
            captured_utc=datetime.now(timezone.utc).isoformat(),
            workers=rows,
        ),
    )
    return destination


def environment(cfg: dict) -> dict[str, str]:
    """Keep each native calculation on its original single-thread CPU environment."""
    amber = Path(cfg["sander"]).parent.parent
    return dict(
        AMBERHOME=str(amber),
        LD_LIBRARY_PATH=str(amber / "lib"),
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        CUDA_VISIBLE_DEVICES="",
    )


def activate(schedule: Path, host: str) -> None:
    """Replace queue control while native Sander finishes its current full block."""
    import signal

    plan, cfg, prepared = context(schedule, host)
    folder = schedule.parent / host
    record = json.loads((folder / "inventory.json").read_text())
    if record["schedule_sha256"] != digest(schedule):
        raise ValueError("Inventory schedule identity changed")
    for row in record["workers"]:
        for role in ("controller", "driver", "native"):
            now = capture(row[role]["pid"])
            if any(
                now[k] != row[role][k] for k in ("start_ticks", "argv", "cwd", "exe")
            ):
                raise ValueError(f"Process changed before handoff: {role}")
    with (folder / "activation_claim").open("x") as stream:
        stream.write(str(os.getpid()) + "\n")
    for row in record["workers"]:
        event(folder, "retire_queue_intent", identity=row["controller"])
        signal_verified(row["controller"], signal.SIGTERM)
        event(folder, "retire_queue_outcome", pid=row["controller"]["pid"])
    for index, row in enumerate(record["workers"]):
        event(
            folder,
            "suspend_window_driver_intent",
            identity=row["driver"],
            native=row["native"],
        )
        suspend_driver(row["native"], row["driver"])
        event(folder, "suspend_window_driver_outcome", pid=row["driver"]["pid"])
        job = folder / f"worker{index + 1:02d}" / "job"
        launch_job(
            f"nyl12_us_coverage_first_{host}_{index + 1:02d}_20261120",
            job,
            prepared.parent,
            [
                cfg["python"],
                plan["script"],
                "worker",
                "--schedule",
                str(schedule),
                "--host",
                host,
                "--index",
                str(index),
            ],
            env_overrides=environment(cfg),
            input_paths=[schedule, folder / "inventory.json"],
        )
        event(folder, "replacement_worker_launched", worker=index, job=str(job))
    write_json(
        folder / "activation.json",
        dict(
            status="workers_launched",
            workers=len(record["workers"]),
            native_blocks_interrupted=0,
            schedule_sha256=digest(schedule),
        ),
    )


def coverage_receipt(root: Path, labels: list[str]) -> dict | None:
    """Require a complete diagnostic-ready first production block in every window."""
    receipts = {}
    required = (
        "acceptance.json",
        "process.json",
        "analysis_request.json",
        "samples_and_cross_bias.npz",
        "observables/observables.csv",
        "observables/observables.json",
    )
    for label in labels:
        stage = root / label / "03_production_block_01_5ps"
        if not all((stage / name).is_file() for name in required):
            return None
        accepted = json.loads((stage / "acceptance.json").read_text())
        process = json.loads((stage / "process.json").read_text())
        if (
            accepted["completed_steps"] != 10000
            or not accepted["finite_energy_temperature"]
        ):
            raise ValueError(f"Unaccepted first block: {label}")
        if process["exit_status"] != 0 or process["ended_utc"] is None:
            raise ValueError(f"Unsuccessful first block process: {label}")
        receipts[label] = digest(stage / "acceptance.json")
    return dict(
        status="all_windows_have_first_block",
        windows=len(labels),
        first_block_acceptance_sha256=receipts,
        released_utc=datetime.now(timezone.utc).isoformat(),
        sampling_converged=False,
    )


def release(schedule: Path, host: str, receipt: dict) -> None:
    """Publish one immutable release after checking the host's own accepted blocks."""
    plan, _, prepared = context(schedule, host)
    labels = [t["label"] for queue in plan["assignments"][host] for t in queue]
    if (
        receipt["windows"] != 63
        or set(receipt["first_block_acceptance_sha256"]) != set(plan["windows"])
        or coverage_receipt(prepared.parent, labels) is None
    ):
        raise ValueError("Cannot release extensions before complete coverage")
    for label in labels:
        path = prepared.parent / label / "03_production_block_01_5ps/acceptance.json"
        if digest(path) != receipt["first_block_acceptance_sha256"][label]:
            raise ValueError(f"First block identity mismatch: {label}")
    gate = schedule.parent / host / "coverage_gate.json"
    if gate.exists():
        raise FileExistsError(gate)
    write_json(gate, receipt)


def coordinate(schedule: Path, host: str, interval: int = 60) -> None:
    """Release both hosts only after all 63 windows publish their first full block."""
    plan, _, prepared = context(schedule, host)
    if host != plan["local_host"]:
        raise ValueError("Global coordinator must run on the analysis host")
    while (receipt := coverage_receipt(prepared.parent, plan["windows"])) is None:
        time.sleep(interval)
    remote = plan["remote_host"]
    remote_argv = [
        plan["host_python"][remote],
        plan["script"],
        "release",
        "--schedule",
        str(schedule),
        "--host",
        remote,
    ]
    subprocess.run(
        [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=15",
            remote,
            shlex.join(remote_argv),
        ],
        input=json.dumps(receipt),
        text=True,
        check=True,
    )
    release(schedule, host, receipt)
