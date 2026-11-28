#!/usr/bin/env python3
"""Prepare or run the isolated four-target relay production continuation."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.relay_targeted_extension import prepare_extension  # noqa: E402
from nyl12_covalent.relay_targeted_extension_runtime import run_extension  # noqa: E402


def main() -> None:
    """Create frozen continuation inputs or execute one reviewed target label."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prepare = sub.add_parser("prepare")
    prepare.add_argument("prepared", type=Path)
    prepare.add_argument("output", type=Path)
    run = sub.add_parser("run")
    run.add_argument("manifest", type=Path)
    run.add_argument("label")
    run.add_argument("--expected-sha256", required=True)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare_extension(args.prepared, args.output)
    else:
        run_extension(args.manifest, args.label, args.expected_sha256)


if __name__ == "__main__":
    main()
