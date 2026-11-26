#!/usr/bin/env python3
"""Render the 300 K relay umbrella pilot without launching Amber jobs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.relay_umbrella_campaign_prep import (  # noqa: E402
    prepare_full_path,
    prepare_pilot,
)

LOGGER = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Prepare immutable per-window Amber inputs and stop before execution."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--campaign",
        type=Path,
        default=ROOT / "qmmm/pa6_flipped_model30/12_nterminal_relay",
    )
    parser.add_argument(
        "--output",
        type=Path,
    )
    parser.add_argument("--full-path", action="store_true")
    args = parser.parse_args(argv)
    output = args.output or args.campaign / (
        "umbrella_sampling/full_path_300K"
        if args.full_path
        else "umbrella_sampling/pilot_300K_2d_v2"
    )
    prepared = (
        prepare_full_path(args.campaign, output)
        if args.full_path
        else prepare_pilot(args.campaign, output)
    )
    LOGGER.info("Prepared relay umbrella inputs: %s", prepared)
    LOGGER.info("No native calculations were launched")
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
