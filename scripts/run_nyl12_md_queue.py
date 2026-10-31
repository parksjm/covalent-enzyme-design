#!/usr/bin/env python3
"""Run one portable, sequential Amber replica queue on one GPU."""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import time
from pathlib import Path

LOGGER = logging.getLogger("run_nyl12_md_queue")


def wait_for_success(status_path: Path, poll_seconds: float) -> None:
    """Wait for a prerequisite scientific job to complete successfully."""
    if poll_seconds <= 0:
        raise ValueError("Poll interval must be positive")
    while True:
        status = json.loads(status_path.read_text(encoding="utf-8"))
        exit_status = status.get("exit_status")
        if exit_status is None:
            time.sleep(poll_seconds)
            continue
        if exit_status != 0:
            raise RuntimeError(f"Prerequisite job failed with status {exit_status}")
        return


def run_queue(campaign: Path, queue_file: Path, gpu: int) -> None:
    """Run every relative replica directory listed in a queue.

    Args:
        campaign: Immutable campaign root containing ``systems``.
        queue_file: Text file with one campaign-relative replica path per line.
        gpu: CUDA device index visible on this host.
    """
    if gpu < 0:
        raise ValueError("GPU index must be nonnegative")
    campaign = campaign.resolve()
    entries = [
        line.strip()
        for line in queue_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not entries:
        raise ValueError(f"Queue is empty: {queue_file}")
    environment = dict(os.environ)
    environment["CUDA_VISIBLE_DEVICES"] = str(gpu)
    for entry in entries:
        replica = (campaign / entry).resolve()
        if campaign not in replica.parents:
            raise ValueError(f"Queue entry escapes campaign root: {entry}")
        runner = replica / "run_replica.sh"
        if not runner.is_file():
            raise FileNotFoundError(runner)
        LOGGER.info("GPU %d starting %s", gpu, entry)
        subprocess.run(
            ["bash", "run_replica.sh"], cwd=replica, env=environment, check=True
        )
        LOGGER.info("GPU %d completed %s", gpu, entry)


def main(argv: list[str] | None = None) -> int:
    """Parse command arguments and run one GPU queue."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--queue", required=True, type=Path)
    parser.add_argument("--gpu", required=True, type=int)
    parser.add_argument("--after-status", type=Path)
    parser.add_argument("--poll-seconds", type=float, default=60.0)
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.after_status is not None:
        LOGGER.info("waiting for successful prerequisite %s", args.after_status)
        wait_for_success(args.after_status, args.poll_seconds)
    run_queue(args.campaign, args.queue, args.gpu)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileNotFoundError,
        RuntimeError,
        subprocess.CalledProcessError,
        ValueError,
    ) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
