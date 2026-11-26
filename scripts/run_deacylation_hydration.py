#!/usr/bin/env python3
"""Prepare or execute a frozen, unbiased deacylation hydration probe."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if (HERE / "nyl12_covalent").is_dir():
    sys.path.insert(0, str(HERE))
else:
    sys.path.insert(0, str(HERE.parent / "src"))

from nyl12_covalent.deacylation_hydration import run  # noqa: E402
from nyl12_covalent.deacylation_hydration_prepare import prepare  # noqa: E402


def main() -> int:
    """Prepare immutable inputs or execute a previously prepared campaign."""
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare-only", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()
    if args.prepare_only:
        if args.audit is None or args.repo is None or args.seed is None:
            parser.error("--prepare-only requires --audit, --repo and --seed")
        prepare(args.campaign, args.audit, args.repo, args.seed)
        return 0
    status = run(args.campaign)
    return 0 if status["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
