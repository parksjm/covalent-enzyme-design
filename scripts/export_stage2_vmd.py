#!/usr/bin/env python3
"""Export one Amber QM/MM optimizer trajectory for stripped VMD inspection."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.stage2_vmd_export import export_stage2  # noqa: E402


def main() -> int:
    """Parse arguments and export the requested stage trajectory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--system-dir", type=Path, required=True)
    parser.add_argument("--inputs-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpptraj", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    export_stage2(
        args.run_dir.resolve(),
        args.system_dir.resolve(),
        args.inputs_json.resolve(),
        args.output.resolve(),
        args.cpptraj.resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
