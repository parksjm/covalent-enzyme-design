#!/usr/bin/env python3
"""Prepare/review/launch bounded native two-CV NEB umbrellas with CPU queues."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.jobs import launch_job  # noqa:E402
from nyl12_covalent.neb_umbrella_campaign import frozen, prepare, write_json  # noqa:E402
from nyl12_covalent.neb_umbrella_runtime import run_window  # noqa:E402


def environment(cfg: dict[str, Any]) -> dict[str, str]:
    """Declare single-thread CPU Sander without changing any running job environment."""
    amber = Path(cfg["sander"]).parent.parent
    return {
        "AMBERHOME": str(amber),
        "LD_LIBRARY_PATH": str(amber / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "CUDA_VISIBLE_DEVICES": "",
    }


def worker(path: Path, index: int, poll: int = 30) -> None:
    """Launch one window at a time; stop this worker on its first diagnostic failure."""
    cfg = frozen(path)
    output = path.parent
    if not 0 <= index < cfg["workers"]:
        raise ValueError("Workerindexoutsidepreparedschedule")
    folder = output / "workers" / f"worker_{index + 1:02d}"
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "claim").open("x") as stream:
        stream.write(str(os.getpid()) + "\n")
    by_label = {r["label"]: r for r in cfg["windows"]}
    tasks = [
        by_label[label] for label in cfg["priority_schedule"][index :: cfg["workers"]]
    ]
    completed = []
    row = None
    try:
        for row in tasks:
            write_json(
                folder / "status.json",
                {"status": "running", "current": row["label"], "completed": completed},
            )
            job = output / row["label"] / "job"
            launch_job(
                f"{cfg['session_prefix']}_{row['label']}",
                job,
                output,
                [
                    cfg["python"],
                    cfg["runtime_script"],
                    "window",
                    "--config",
                    str(path),
                    "--window",
                    row["label"],
                ],
                env_overrides=environment(cfg),
                input_paths=[path, *[Path(p) for p in cfg["prepared_input_sha256"]]],
            )
            while True:
                raw = (job / "status.json").read_bytes()
                status = json.loads(raw) if raw.endswith(b"\n") else None
                if status is not None and status.get("ended_utc") is not None:
                    break
                time.sleep(poll)
            if status["exit_status"] != 0:
                raise RuntimeError(
                    f"Window{row['label']}failed; inspectitsretainedstatus/logs"
                )
            receipt = json.loads(
                (output / row["label"] / "window_status.json").read_text()
            )
            if receipt["status"] != "completed":
                raise ValueError("Missingcompletewindowreceipt")
            completed.append(row["label"])
        write_json(
            folder / "status.json",
            {
                "status": "completed",
                "completed": completed,
                "sampling_convergence": "Notasserted",
            },
        )
    except BaseException as error:
        write_json(
            folder / "status.json",
            {
                "status": "failed",
                "current": row["label"] if row else None,
                "completed": completed,
                "error": f"{type(error).__name__}: {error}",
            },
        )
        raise


def launch(path: Path) -> None:
    """Start at most eight independent CPU workers after explicit root review."""
    cfg = frozen(path)
    output = path.parent
    for index in range(cfg["workers"]):
        folder = output / "workers" / f"worker_{index + 1:02d}"
        folder.mkdir(parents=True, exist_ok=True)
        launch_job(
            f"{cfg['session_prefix']}_worker{index + 1:02d}",
            folder / "job",
            output,
            [
                cfg["python"],
                cfg["runtime_script"],
                "worker",
                "--config",
                str(path),
                "--worker",
                str(index),
            ],
            env_overrides=environment(cfg),
            input_paths=[path, *[Path(p) for p in cfg["prepared_input_sha256"]]],
        )
    write_json(
        output / "status.json",
        {
            "status": "workers_launched",
            "workers": cfg["workers"],
            "window_count": len(cfg["windows"]),
            "production_ps_minimum": 30,
        },
    )


def main() -> None:
    """Require explicit prepare/launch actions; window is the internal per-job action."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "launch", "worker", "window"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--worker", type=int)
    parser.add_argument("--window")
    args = parser.parse_args()
    path = args.config.resolve(strict=True)
    if args.action == "prepare":
        prepared = prepare(path)
        logging.info("Prepared campaign: %s", prepared)
    elif args.action == "launch":
        launch(path)
    elif args.action == "worker":
        if args.worker is None:
            parser.error("worker requires --worker index")
        worker(path, args.worker)
    else:
        if args.window is None:
            parser.error("window requires --window label")
        run_window(path, args.window)


if __name__ == "__main__":
    main()
