#!/usr/bin/env python3
"""Export a portable stripped VMD movie of deacylation hydration sampling."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.hydration_vmd_exporter import export  # noqa: E402


def main() -> int:
    """Parse arguments and export the hydration visualization."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--water-summary", type=Path, required=True)
    parser.add_argument("--cpptraj", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    result = export(
        args.source.resolve(),
        args.output.resolve(),
        args.cpptraj.resolve(),
        args.water_summary.resolve(),
    )
    logging.getLogger(__name__).info("Export result: %s", result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
