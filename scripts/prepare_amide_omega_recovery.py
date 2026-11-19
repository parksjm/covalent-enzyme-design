#!/usr/bin/env python3
"""Prepare non-overwriting recovery inputs for the PA6 amide omega scan."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.amide_omega_recovery import (  # noqa: E402
    prepare_amide_omega_recovery,
)

LOGGER = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Prepare recovery Gaussian inputs from validated original scan files.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-scan-dir", required=True, type=Path)
    parser.add_argument("--source-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    prepare_amide_omega_recovery(
        args.source_scan_dir, args.source_manifest, args.output_dir
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
