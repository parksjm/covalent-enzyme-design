#!/usr/bin/env python3
"""Wait for twelve matched WT MD replicas and execute a frozen analysis workflow."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.matched_wt_analysis_queue import run_queue  # noqa: E402


def main() -> None:
    """Run the configured completion gate in a new status directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    run_queue(json.loads(args.config.read_text()), args.output)


if __name__ == "__main__":
    main()
