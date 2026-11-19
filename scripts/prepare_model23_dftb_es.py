#!/usr/bin/env python3
"""Prepare the fully unrestrained Model 23 DFTB3/MM E-S optimization."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.es_dftb_unrestrained import (  # noqa: E402
    prepare_dftb_optimization,
)


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-config", type=Path, required=True)
    parser.add_argument("--sander", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--active-mask", default="@*")
    parser.add_argument("--case-label", default="Model 23")
    return parser.parse_args()


def main() -> int:
    """Prepare and report the audited DFTB3/MM command."""
    args = parse_args()
    metadata = prepare_dftb_optimization(
        args.region_config,
        args.sander,
        args.output,
        args.active_mask,
        args.case_label,
    )
    LOGGER.info("Prepared command: %s", " ".join(metadata["command"]))
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
