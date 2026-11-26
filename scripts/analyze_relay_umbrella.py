#!/usr/bin/env python3
"""Analyze immutable compact accepted blocks from both relay umbrella charts."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.relay_umbrella_analysis_io import analyze_request
from nyl12_covalent.relay_umbrella_analysis_reports import compare_reports, plot_report


def main() -> None:
    """Run one analysis request; fail explicitly if source or input checks fail."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--previous", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Retain prior reports; choose a new output path")
    result = analyze_request(args.request.resolve(strict=True))
    if args.previous:
        result["gates"]["block_stability"] = compare_reports(
            result, json.loads(args.previous.read_text())
        )
    if args.plot:
        plot_report(result, args.plot)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    logging.info("Analysis status: %s; report: %s", result["status"], args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
