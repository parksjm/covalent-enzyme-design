"""Summarize cpptraj distance tables from a proton-state audit."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.proton_state_audit import write_summary  # noqa: E402


def main() -> None:
    """Parse command-line arguments and write a JSON stage summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("tables", nargs="+", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    write_summary(args.tables, args.output)


if __name__ == "__main__":
    main()
