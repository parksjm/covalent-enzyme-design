#!/usr/bin/env python3
"""Prepare bidirectional relaxed PA6 Gaussian torsion/improper scans."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_torsion_scans import prepare_relaxed_scans  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Prepare all declared scan jobs.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    prepare_relaxed_scans(args.geometry_dir, args.output_dir)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
