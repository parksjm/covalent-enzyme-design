"""Validate or execute a frozen matched six-case Violet Amber plan."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.matched_original_amber import run_campaign, validate_plan


def main() -> None:
    """Parse an explicit plan and optional single-case worker selection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--case")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    if args.validate_only:
        validate_plan(plan)
        logging.info("Matched plan input/source pins accepted; no calculations run")
        return
    if args.output is None:
        parser.error("--output is required for execution")
    report = run_campaign(plan, args.output, case_name=args.case)
    if report["failure_count"]:
        raise SystemExit(1)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
