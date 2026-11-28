"""Wait for all original relay states at60ps, then execute six reviewed checks."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import logging
from pathlib import Path
import subprocess
import time
from typing import Any

from .neb_umbrella_campaign import digest, write_json
from .relay_mixing_checks import TARGETS
from .relay_targeted import verify_target

LOGGER = logging.getLogger(__name__)
MAIN_LABELS = tuple(
    [f"activation_full_{i:03d}" for i in range(1, 32)]
    + [f"delivery_full_{i:03d}" for i in range(1, 27)]
)


def main_gate(base: Path, active_main_workers: list[int]) -> dict[str, Any]:
    """Return readiness only for57 accepted original stage14 native checkpoints."""
    pending = []
    for label in MAIN_LABELS:
        stage = base / label / "14_production_5ps"
        paths = [
            stage / n for n in ("acceptance.json", "process.json", "restart.ncrst")
        ]
        if not all(p.is_file() for p in paths):
            pending.append(label)
            continue
        acceptance = json.loads(paths[0].read_text())
        process = json.loads(paths[1].read_text())
        if not (
            acceptance.get("accepted") is True
            and acceptance.get("completed_steps") == 10000
            and acceptance.get("every_integration_step_SCC_converged") is True
            and acceptance.get("frozen_outer_max_displacement_A") == 0
            and process.get("exit_status") == 0
            and process.get("ended_utc")
            and acceptance.get("checkpoint_sha256") == digest(paths[2])
        ):
            pending.append(label)
    return dict(
        ready=not pending and not active_main_workers,
        original_main_count=57,
        accepted60ps_count=57 - len(pending),
        pending_labels=pending,
        active_main_worker_pids=active_main_workers,
    )


def main_worker_pids(base: Path) -> list[int]:
    """Inspect live processes by direct main-window working directory and argv."""
    found = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            argv = (proc / "cmdline").read_bytes().replace(b"\0", b" ").decode()
            cwd = (proc / "cwd").resolve(strict=True)
        except (PermissionError, FileNotFoundError, ProcessLookupError):
            continue
        if (
            "sander" in argv
            and cwd.parent.parent == base
            and cwd.parent.name in MAIN_LABELS
        ):
            found.append(int(proc.name))
        elif "run_relay_umbrella" in argv and str(base) in argv:
            found.append(int(proc.name))
    return found


def validate_method(cfg: dict[str, Any]) -> None:
    """Reject changed charge, masks, temperature, model size or production horizon."""
    system = cfg["system"]
    expected = dict(
        full_atoms=104922,
        active_atoms=4165,
        qm_real_atoms=98,
        qm_mm_link_atoms=8,
        qm_charge=-2,
        qm_method="DFTB3/3ob-3-1",
    )
    if any(system.get(k) != v for k, v in expected.items()):
        raise ValueError("Reviewed method/system inventory changed")
    if cfg["temperature_K"] != 300 or cfg["production_horizon_ps"] != 30:
        raise ValueError("Reviewed thermal protocol changed")


def execute_task(task: dict[str, Any], output: Path) -> None:
    """Run one task with its own process status and log, retaining failures."""
    import os

    label = task["label"]
    status = output / f"{label}.json"
    write_json(status, dict(label=label, status="running", argv=task["argv"]))
    with (output / f"{label}.log").open("xb") as log:
        result = subprocess.run(
            task["argv"],
            cwd=task["cwd"],
            env={**os.environ, **task["environment"]},
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    write_json(
        status,
        dict(
            label=label,
            status="completed" if result.returncode == 0 else "failed",
            exit_status=result.returncode,
        ),
    )
    if result.returncode:
        raise RuntimeError(f"Independent check {label} failed: {result.returncode}")


def run_queue(tasks_path: Path, output: Path, poll_seconds: float = 300) -> None:
    """Persistently gate six unique tasks, then launch exactly six CPU workers."""
    data = json.loads(tasks_path.read_text())
    tasks = data["tasks"]
    if len(tasks) != 6 or {t["label"] for t in tasks} != set(TARGETS):
        raise ValueError("Queue requires exactly six unique approved labels")
    prepared = Path(data["prepared"])
    cfg = verify_target(prepared, data["prepared_sha256"])
    validate_method(cfg)
    output.mkdir(exist_ok=False)
    for task in tasks:
        write_json(
            output / f"{task['label']}.json", dict(label=task["label"], status="queued")
        )
    base = Path(cfg["parent_prepared"]).parent
    while True:
        gate = main_gate(base, main_worker_pids(base))
        write_json(output / "gate.json", gate)
        LOGGER.info(
            "Main accepted60ps=%d/57 workers=%d",
            gate["accepted60ps_count"],
            len(gate["active_main_worker_pids"]),
        )
        if gate["ready"]:
            break
        time.sleep(poll_seconds)
    verify_target(prepared, data["prepared_sha256"])
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = [pool.submit(execute_task, task, output) for task in tasks]
        for future in futures:
            future.result()
    write_json(
        output / "queue_completion.json",
        dict(
            process_completed=True,
            scientific_acceptance="See individual native window_status/acceptance receipts",
        ),
    )


def main() -> None:
    """Run the explicit host-specific six-check deferred queue."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tasks", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run_queue(args.tasks, args.output)


if __name__ == "__main__":
    main()
