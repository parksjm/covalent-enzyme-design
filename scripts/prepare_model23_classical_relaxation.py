#!/usr/bin/env python3
"""Prepare a fully unrestrained classical relaxation of PA6 Model 23."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.es_classical_relaxation import prepare_relaxation  # noqa: E402


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--engine", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case-label", default="Model 23")
    return parser.parse_args()


def main() -> int:
    """Write the relaxation input and provenance record."""
    args = parse_args()
    metadata = prepare_relaxation(
        args.topology,
        args.coordinates,
        args.engine,
        args.output,
        args.case_label,
    )
    LOGGER.info("Prepared command: %s", " ".join(metadata["command"]))
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
