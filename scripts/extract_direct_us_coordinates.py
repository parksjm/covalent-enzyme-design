#!/usr/bin/env python3
"""Augment an explicit host-local list of accepted direct-US production stages."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import logging
from pathlib import Path

from nyl12_covalent.direct_us_coordinates import augment_stage


def _extract(arguments: tuple[str, str]) -> dict:
    return augment_stage(Path(arguments[0]), Path(arguments[1]))


def main() -> None:
    """Run only planned stage extractions and retain a complete process summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--workers", required=True, type=int)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8 or args.summary.exists():
        raise ValueError("Require1–8workers and a new summary path")
    plan = json.loads(args.plan.read_text())
    stages = plan["stages"]
    if not stages or len({row["stage"] for row in stages}) != len(stages):
        raise ValueError("Plan must explicitly name unique host-local stages")
    arguments = [(row["stage"], row["output"]) for row in stages]
    if args.workers == 1:
        results = [_extract(item) for item in arguments]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(_extract, arguments))
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(
        json.dumps(
            dict(
                host=plan["host"],
                stages=results,
                count=len(results),
                scientific_jobs_changed=False,
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    logging.info(
        "Completed %s coordinate augmentations on %s", len(results), plan["host"]
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
