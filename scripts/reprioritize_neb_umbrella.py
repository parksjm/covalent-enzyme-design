#!/usr/bin/env python3
"""Run an explicitly approved, immutable coverage-first umbrella schedule."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def main() -> None:
    """Require exact scheduling identity, execution host, and an explicit action."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=("inventory", "activate", "worker", "coordinate", "release")
    )
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--host", required=True)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    from nyl12_covalent import neb_umbrella_coverage as coverage

    if args.action == "worker":
        from nyl12_covalent.neb_umbrella_coverage_worker import worker

        if args.index is None:
            parser.error("worker requires --index")
        worker(args.schedule, args.host, args.index)
    elif args.action == "release":
        coverage.release(args.schedule, args.host, json.load(sys.stdin))
    else:
        getattr(coverage, args.action)(args.schedule, args.host)


if __name__ == "__main__":
    main()
