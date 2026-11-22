"""Command-line entry point for native NEB-US movie-part export."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.neb_us_movie import export_part


def main() -> None:
    """Parse arguments and export selected umbrella windows."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--windows", nargs="+", required=True)
    parser.add_argument("--cpptraj", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frames-per-block", type=int, default=5)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    export_part(
        args.campaign, args.windows, args.cpptraj, args.output, args.frames_per_block
    )


if __name__ == "__main__":
    main()
