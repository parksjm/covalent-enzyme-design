#!/usr/bin/env python3
"""Summarize completed matched WT relay/native geometry records."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.matched_wt_relays import summarize  # noqa: E402


def main() -> None:
    """Parse explicit completed analysis inputs and write a new joined summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--occupancy-dir", type=Path, required=True)
    parser.add_argument("--screen-dir", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    summarize(args.config, args.occupancy_dir, args.screen_dir, args.output)
    logging.info("Wrote matched WT joined summary: %s", args.output)


if __name__ == "__main__":
    main()
