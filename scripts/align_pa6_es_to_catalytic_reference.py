"""Align converged PA6 E-S coordinates to a CIF catalytic scaffold."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.catalytic_alignment import align_catalytic_site  # noqa: E402


def main() -> None:
    """Run the catalytic-site alignment from explicit input and output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_cif", type=Path)
    parser.add_argument("topology", type=Path)
    parser.add_argument("starting_restart", type=Path)
    parser.add_argument("restart", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)
    align_catalytic_site(
        args.reference_cif,
        args.topology,
        args.starting_restart,
        args.restart,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
