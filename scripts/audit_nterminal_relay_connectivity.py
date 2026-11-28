#!/usr/bin/env python3
"""Run or watch the read-only N-terminal relay connectivity audit."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.relay_connectivity_report import refresh_watch


def main() -> None:
    """Parse CLI options and run one audit or a bounded refresh watcher."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval-seconds", type=float, default=30.0)
    parser.add_argument("--max-seconds", type=float, default=43200.0)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    output = args.output or args.campaign / "connectivity_audit"
    refresh_watch(
        args.campaign,
        output,
        interval_s=args.interval_seconds,
        max_seconds=args.max_seconds,
        once=not args.watch,
    )


if __name__ == "__main__":
    main()
