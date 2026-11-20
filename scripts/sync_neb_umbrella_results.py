#!/usr/bin/env python3
"""Freeze or run compact remote completed-block synchronization without raw MD."""

from __future__ import annotations

import argparse
import fcntl
import logging
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_umbrella_sync_transport import (  # noqa: E402
    freeze_sources,
    sync_once,
    verify_sources,
)


def main() -> None:
    """Require explicit remote ownership, frozen source and one exclusive watcher."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--remote", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--freeze", action="store_true")
    mode.add_argument("--once", action="store_true")
    mode.add_argument("--interval", type=int)
    parser.add_argument("--duration-hours", type=float, default=168)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    path = args.plan.resolve(strict=True)
    if args.freeze:
        logging.info("Frozen synchronizer: %s", freeze_sources(path, args.remote))
        return
    if args.interval is not None and (args.interval < 1 or args.duration_hours < 72):
        parser.error(
            "watch requires positive interval and duration of at least 72 hours"
        )
    with (path.parent / f"sync_{args.remote}.lock").open("a") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        deadline = time.monotonic() + args.duration_hours * 3600
        while True:
            verify_sources(Path(__file__))
            counts = sync_once(path, args.remote)
            logging.info("Synced completed compact artifacts: %s", counts)
            if args.once or time.monotonic() >= deadline:
                break
            time.sleep(args.interval)


if __name__ == "__main__":
    main()
