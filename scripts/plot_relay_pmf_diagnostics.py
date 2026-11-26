#!/usr/bin/env python3
"""Render standardized relay PMF uncertainty and overlap diagnostics."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.relay_pmf_plotting import render_diagnostics


def main() -> None:
    """Run the report-only plotting command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--uncertainty", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stem", required=True)
    parser.add_argument(
        "--interval-method",
        choices=("analytic", "block-bootstrap"),
        required=True,
        help="Select the interval shown; there is no implicit fallback between methods.",
    )
    args = parser.parse_args()
    render_diagnostics(
        args.report,
        args.uncertainty,
        args.output_dir,
        args.stem,
        args.interval_method,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
