#!/usr/bin/env python3
"""Prepare and run the reviewed direct-US 30-to-60 ps continuation."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.direct_us_continuation import prepare_continuation  # noqa: E402
from nyl12_covalent.direct_us_continuation_scheduler import run_worker  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser(
        "prepare", help="pin plan and continuation inputs; do not launch"
    )
    prepare.add_argument("--campaign", type=Path, required=True)
    prepare.add_argument("--control", type=Path, required=True)
    prepare.add_argument("--workers", type=int, default=40)
    worker = commands.add_parser("worker", help="run dynamic 5 ps block-boundary work")
    worker.add_argument("--plan", type=Path, required=True)
    worker.add_argument("--worker-id", type=int, required=True)
    worker.add_argument("--idle-seconds", type=float, default=2.0)
    return parser


def main() -> int:
    """Parse one explicit operation and emit concise operational logs."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parser().parse_args()
    if args.command == "prepare":
        plan = prepare_continuation(args.campaign, args.control, args.workers)
        LOGGER.info("Prepared immutable continuation plan: %s", plan)
        return 0
    result = run_worker(args.plan, args.worker_id, args.idle_seconds)
    LOGGER.info("Worker %d finished: %s", args.worker_id, result["status"])
    return 0 if result["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
