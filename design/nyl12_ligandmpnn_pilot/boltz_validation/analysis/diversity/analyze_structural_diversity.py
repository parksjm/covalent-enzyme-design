#!/usr/bin/env python3
"""Run the retained offline Nyl12 structural-diversity analysis."""

from __future__ import annotations

import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT))

from diversity_reporting import run_analysis  # noqa: E402


def main() -> int:
    """Analyze retained raw/relaxed structures and write compact results."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    directory = Path(__file__).resolve().parent
    root = directory.parents[4]
    boltz_root = root / "design/nyl12_ligandmpnn_pilot/boltz_validation"
    run_analysis(root, boltz_root, directory, Path(__file__))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
