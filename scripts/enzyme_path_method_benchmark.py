#!/usr/bin/env python3
"""Prepare, convert, or analyze the enzyme-path method benchmark."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.enzyme_path_method_analysis import (  # noqa: E402
    analyze,
    correct_periodic_profile,
)
from nyl12_covalent.enzyme_path_method_benchmark import (  # noqa: E402
    prepare_benchmark,
    write_xtb_embedding,
)


def main() -> None:
    """Dispatch an explicit benchmark operation."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--coordinate-drive", required=True, type=Path)
    prepare.add_argument("--topology", required=True, type=Path)
    prepare.add_argument("--output", required=True, type=Path)
    convert = subparsers.add_parser("convert")
    convert.add_argument("--gaussian-input", required=True, type=Path)
    convert.add_argument("--output", required=True, type=Path)
    analysis = subparsers.add_parser("analyze")
    analysis.add_argument("--input-root", required=True, type=Path)
    analysis.add_argument("--output", required=True, type=Path)
    correction = subparsers.add_parser("correct-periodic")
    correction.add_argument("--matched-energies", required=True, type=Path)
    correction.add_argument("--coordinate-drive", required=True, type=Path)
    correction.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    if args.command == "prepare":
        prepare_benchmark(args.coordinate_drive, args.topology, args.output)
    elif args.command == "convert":
        write_xtb_embedding(args.gaussian_input, args.output)
    elif args.command == "analyze":
        analyze(args.input_root, args.output)
    else:
        correct_periodic_profile(
            args.matched_energies, args.coordinate_drive, args.output
        )


if __name__ == "__main__":
    main()
