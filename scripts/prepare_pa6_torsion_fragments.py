#!/usr/bin/env python3
"""Prepare PA6 torsion-benchmark fragment seeds and Gaussian inputs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_torsion_fragments import prepare_torsion_fragments  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Prepare declared benchmark fragments.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    prepare_torsion_fragments(args.output_dir)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, RuntimeError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
