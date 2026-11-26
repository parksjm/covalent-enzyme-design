#!/usr/bin/env python3
"""Prepare or execute the sequential hydration MD continuation."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if (HERE.parent / "nyl12_covalent").is_dir():
    sys.path.insert(0, str(HERE.parent))
else:
    sys.path.insert(0, str(HERE.parent / "src"))

from nyl12_covalent.deacylation_hydration_extension import (  # noqa: E402
    prepare_extension,
    run_extension,
)


def main() -> int:
    """Prepare immutable inputs or run all four reviewed segments."""
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare-only", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--repo", type=Path)
    args = parser.parse_args()
    if args.prepare_only:
        if args.pilot is None or args.repo is None:
            parser.error("--prepare-only requires --pilot and --repo")
        prepare_extension(args.campaign, args.pilot, args.repo)
        return 0
    result = run_extension(args.campaign)
    return 0 if result["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
