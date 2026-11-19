#!/usr/bin/env python3
"""Wait for one GPU to become idle, then run one declared command."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
import subprocess
import time

LOGGER = logging.getLogger(__name__)


def parse_pids(output: str) -> list[int]:
    """Parse the no-header PID output from ``nvidia-smi``."""
    pids = []
    for line in output.splitlines():
        value = line.strip()
        if value and value not in {"N/A", "[Not Supported]"}:
            pids.append(int(value))
    return pids


def active_pids(nvidia_smi: Path, gpu: int) -> list[int]:
    """Return compute-process PIDs currently assigned to one GPU."""
    result = subprocess.run(
        [
            str(nvidia_smi),
            "-i",
            str(gpu),
            "--query-compute-apps=pid",
            "--format=csv,noheader,nounits",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return parse_pids(result.stdout)


def wait_and_run(
    nvidia_smi: Path,
    gpu: int,
    poll_seconds: int,
    command: list[str],
) -> int:
    """Wait without preempting active work, then run the command on the GPU."""
    if gpu < 0 or poll_seconds < 1 or not command:
        raise ValueError(
            "Nonnegative GPU, positive poll interval, and command required"
        )
    polls = 0
    while True:
        pids = active_pids(nvidia_smi, gpu)
        if not pids:
            break
        if polls == 0 or polls % 10 == 0:
            LOGGER.info("GPU %d occupied by PIDs %s; waiting", gpu, pids)
        polls += 1
        time.sleep(poll_seconds)
    environment = dict(os.environ)
    environment["CUDA_VISIBLE_DEVICES"] = str(gpu)
    LOGGER.info("GPU %d is idle; launching: %s", gpu, command)
    return subprocess.run(command, env=environment, check=False).returncode


def main() -> None:
    """Parse CLI arguments, wait for an idle GPU, and propagate exit status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nvidia-smi", type=Path, default=Path("/usr/bin/nvidia-smi"))
    parser.add_argument("--gpu", required=True, type=int)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    raise SystemExit(
        wait_and_run(args.nvidia_smi, args.gpu, args.poll_seconds, command)
    )


if __name__ == "__main__":
    main()
