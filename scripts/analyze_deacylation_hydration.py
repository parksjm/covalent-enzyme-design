#!/usr/bin/env python3
"""Run descriptive water-geometry analysis for saved deacylation trajectories."""

from __future__ import annotations
import argparse
import logging
from pathlib import Path
from nyl12_covalent.deacylation_water_analysis import analyze


def main() -> None:
    """Parse paths and execute the descriptive analysis."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base",
        type=Path,
        default=Path(
            "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
        ),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    analyze(args.output or args.base / "analysis", args.base)


if __name__ == "__main__":
    main()
