#!/usr/bin/env python3
"""Compare matched relay PMFs from blockwise and continuous-window selections."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.relay_pmf_comparison import compare_reports


def main() -> None:
    """Run report-only comparison plotting."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--continuous", type=Path, required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument(
        "--interval-method", choices=("analytic", "bootstrap"), required=True
    )
    args = parser.parse_args()
    compare_reports(
        args.baseline, args.continuous, args.output_stem, args.interval_method
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
