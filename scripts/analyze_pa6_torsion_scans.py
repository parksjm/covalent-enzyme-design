#!/usr/bin/env python3
"""Analyze bidirectional Gaussian PA6 relaxed torsion scans."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_torsion_scan_analysis import analyze_manifest  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Analyze a Gaussian relaxed-scan manifest.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    report = analyze_manifest(args.manifest, args.output_dir)
    LOGGER.info(
        "Wrote analysis for %d scans to %s", len(report["scans"]), args.output_dir
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
