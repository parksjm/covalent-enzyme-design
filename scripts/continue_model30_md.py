#!/usr/bin/env python3
"""Queue a 90 ns continuation after an accepted 10 ns model30 pilot."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import re
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.run_model30_md_replica import (  # noqa: E402
    digest,
    stage_acceptance,
    write_json,
)

LOGGER = logging.getLogger(__name__)


def continuation_input(text: str, seed: int) -> str:
    """Extend the exact unrestrained 10 ns protocol to a new 90 ns segment.

    Args:
        text: Existing 10 ns production Amber input.
        seed: Positive distinct Langevin seed for the continuation.

    Returns:
        Same protocol with 45 million steps and the explicit new seed.
    """
    if seed <= 0:
        raise ValueError("Continuation seed must be positive")
    for key, expected in (
        ("nstlim", "5000000"),
        ("dt", "0.002"),
        ("irest", "1"),
        ("ntx", "5"),
        ("ntr", "0"),
        ("nmropt", "0"),
    ):
        match = re.search(rf"\b{key}\s*=\s*([^,\s]+)", text)
        if not match or match[1] != expected:
            raise ValueError(f"Unexpected source production control {key}")
    text, steps = re.subn(r"\bnstlim\s*=\s*5000000\b", "nstlim=45000000", text)
    text, seeds = re.subn(r"\big\s*=\s*\d+", f"ig={seed}", text)
    if steps != 1 or seeds != 1:
        raise ValueError("Ambiguous production steps or seed")
    return "06_production_90ns checkpoint continuation\n" + text.split("\n", 1)[1]


def main() -> None:
    """Wait for accepted 10 ns, then append 90 ns without changing old outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("replica-dir", "source-job-dir", "engine"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    folder = args.replica_dir
    source_status = args.source_job_dir / "status.json"
    if not folder.is_dir() or not source_status.is_file():
        raise FileNotFoundError("Replica directory or source job status is missing")
    state_path = folder / "extension_100ns.json"
    if state_path.exists():
        raise FileExistsError(state_path)
    record = {
        "status": "queued_waiting_for_accepted_10ns",
        "additional_ns": 90,
        "total_production_ns": 100,
        "seed": args.seed,
        "queued_utc": datetime.now(timezone.utc).isoformat(),
        "source_job_status": str(source_status),
        "source_sha256": {
            str(p): digest(p)
            for p in (
                Path(__file__),
                Path(__file__).parent / "run_model30_md_replica.py",
                args.engine,
            )
        },
    }
    write_json(state_path, record)
    LOGGER.info("Queued90ns extension; waiting for accepted10ns in %s", folder)
    while not (folder / "completion.json").is_file():
        status = json.loads(source_status.read_text())
        if status["exit_status"] is not None:
            raise RuntimeError(
                f"Source10ns job ended without accepted completion: {status['exit_status']}"
            )
        time.sleep(30)
    source_output = folder / "05_production_10ns.out"
    source_acceptance = stage_acceptance(source_output, 5_000_000)
    checkpoint = folder / "05_production_10ns.rst7"
    topology = folder / "system.prmtop"
    source_input = folder / "05_production_10ns.in"
    mdin = folder / "06_production_90ns.in"
    if mdin.exists():
        raise FileExistsError(mdin)
    mdin.write_text(continuation_input(source_input.read_text(), args.seed))
    stem = folder / "06_production_90ns"
    command = [
        str(args.engine),
        "-O",
        "-i",
        str(mdin),
        "-o",
        str(stem) + ".out",
        "-p",
        str(topology),
        "-c",
        str(checkpoint),
        "-r",
        str(stem) + ".rst7",
        "-x",
        str(stem) + ".nc",
        "-inf",
        str(stem) + ".mdinfo",
    ]
    record.update(
        status="running90ns_checkpoint_continuation",
        command=command,
        source_checkpoint_acceptance=source_acceptance,
        started_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256={
            str(p): digest(p) for p in (topology, checkpoint, mdin, source_input)
        },
    )
    write_json(state_path, record)
    LOGGER.info("Accepted10ns; continuing90ns from %s", checkpoint)
    with Path(str(stem) + ".engine.log").open("w") as handle:
        result = subprocess.run(
            command, cwd=folder, stdout=handle, stderr=subprocess.STDOUT, check=False
        )
    record.update(
        returncode=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    if result.returncode:
        record["status"] = "engine_failed"
        write_json(state_path, record)
        raise RuntimeError(f"90ns extension failed: {result.returncode}")
    acceptance = stage_acceptance(Path(str(stem) + ".out"), 45_000_000)
    acceptance.update(
        total_production_ns=100, restart_sha256=digest(Path(str(stem) + ".rst7"))
    )
    write_json(Path(str(stem) + ".acceptance.json"), acceptance)
    record["status"] = "100ns_complete_requires_projected_phi_native_sampling_analysis"
    write_json(state_path, record)
    LOGGER.info("100ns total production completed")


if __name__ == "__main__":
    main()
