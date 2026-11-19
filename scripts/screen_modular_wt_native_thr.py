#!/usr/bin/env python3
"""Rank joint-eligible WT MD frames by native Thr233 displacement."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.wt_native_thr_screen import screen  # noqa: E402


def main() -> None:
    """Parse arguments and screen configured WT MD frames."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--top-k", type=int, required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(message)s",
    )
    screen(args.config, args.output, args.cpptraj, args.top_k)


if __name__ == "__main__":
    main()
