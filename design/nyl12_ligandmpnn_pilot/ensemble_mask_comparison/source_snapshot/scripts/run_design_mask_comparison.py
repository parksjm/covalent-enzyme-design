#!/usr/bin/env python3
"""Run a manifest-frozen LigandMPNN design-mask comparison."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.nyl12_covalent.design_mask_generation import run_comparison  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_comparison(args.manifest, args.output, args.workers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
