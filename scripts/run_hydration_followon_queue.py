#!/usr/bin/env python3
"""Run only after lead review through the persistent job launcher."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nyl12_covalent.hydration_followon_queue import run  # noqa: E402


def main() -> None:
    """Parse the immutable reviewed plan and run its bounded queue."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run(args.plan)


if __name__ == "__main__":
    main()
