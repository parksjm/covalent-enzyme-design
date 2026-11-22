#!/usr/bin/env python3
"""Render analysis JSON using an explicitly selected existing Matplotlib runtime."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_umbrella_analysis_plots import export  # noqa:E402


def main() -> None:
    """Render one immutable JSON generation without any scientific calculation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", required=True, type=Path)
    args = parser.parse_args()
    export(args.analysis.parent, json.loads(args.analysis.read_text()))


if __name__ == "__main__":
    main()
