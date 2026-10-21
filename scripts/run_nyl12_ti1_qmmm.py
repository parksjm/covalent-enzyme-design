#!/usr/bin/env python3
"""Run a Nyl12 TI1 QM/MM attempt using explicit scientific inputs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.ti1_attempt import run_attempt  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    """Build the command-line parser for a TI1 QM/MM attempt."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--amber-bin", required=True, type=Path)
    parser.add_argument("--target", required=True, choices=("pa6_l2", "poly_xama_maxa"))
    parser.add_argument("--max-cycles", type=int, default=300)
    parser.add_argument("--checkpoint-attempt", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Parse explicit inputs and run the selected TI1 QM/MM attempt.

    Args:
        argv: Optional argument vector; defaults to process arguments.

    Returns:
        Process exit status.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    result = run_attempt(
        topology=args.topology,
        coordinates=args.coordinates,
        reference=args.reference,
        output=args.output,
        amber_bin=args.amber_bin,
        target=args.target,
        max_cycles=args.max_cycles,
        checkpoint_attempt=args.checkpoint_attempt,
    )
    LOGGER.info("TI1 QM/MM attempt completed: %s", result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
