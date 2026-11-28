#!/usr/bin/env python3
"""Prepare targeted relay inputs or execute one caller-pinned target label."""

import argparse
from pathlib import Path

from nyl12_covalent.relay_targeted import prepare_targets
from nyl12_covalent.relay_targeted_runtime import run_target


def main() -> None:
    """Dispatch explicit preparation or one bounded native target execution."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--plan", required=True, type=Path)
    run = commands.add_parser("run")
    run.add_argument("--prepared", required=True, type=Path)
    run.add_argument("--label", required=True)
    run.add_argument("--sha256", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare_targets(args.plan)
    else:
        run_target(args.prepared, args.label, args.sha256)


if __name__ == "__main__":
    main()
