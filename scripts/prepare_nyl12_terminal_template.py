"""Prepare the neutral threonine methylamide charge-template inputs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.terminal_template import (  # noqa: E402
    prepare_templates,
)

LOGGER = logging.getLogger(__name__)


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments for template preparation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new output directory")
    parser.add_argument(
        "--geometry-method",
        required=True,
        choices=("wb97xd-def2tzvp", "b3lyp-d3bj-def2tzvp"),
        help="explicit method for conformer geometry optimization",
    )
    parser.add_argument(
        "--seed",
        dest="seeds",
        type=int,
        action="append",
        help="ETKDGv3 seed; repeat once for each conformer (default: two reviewed seeds)",
    )
    return parser.parse_args()


def main() -> None:
    """Prepare conformer inputs and report the created output directory."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    arguments = _parse_args()
    seeds = tuple(arguments.seeds) if arguments.seeds else (20261020, 20261021)
    output = prepare_templates(
        arguments.output,
        geometry_method=arguments.geometry_method,
        seeds=seeds,
    )
    LOGGER.info("Prepared %s", output)


if __name__ == "__main__":
    main()
