"""Prepare coordinate-matched WT protonation states; never launch MD."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.wt_matched_preparation import prepare


def main() -> None:
    """Run explicit matched preparation into a new directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    prepare(args.root.resolve(strict=True), args.output.resolve())


if __name__ == "__main__":
    main()
