#!/usr/bin/env python3
"""Launch a logged command in a detached tmux session."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.jobs import launch_job  # noqa: E402

LOGGER = logging.getLogger("run_nyl12_job")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    launch = subparsers.add_parser("launch", help="launch a new tmux job")
    launch.add_argument("--session", required=True)
    launch.add_argument("--output-dir", required=True, type=Path)
    launch.add_argument("--cwd", required=True, type=Path)
    launch.add_argument("--env", action="append", default=[], metavar="KEY=VALUE")
    launch.add_argument("--input", action="append", default=[], type=Path)
    launch.add_argument("command", nargs=argparse.REMAINDER)
    status = subparsers.add_parser("status", help="show a job status.json")
    status.add_argument("output_dir", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the launch or status command.

    Args:
        argv: Optional argument vector; defaults to process arguments.

    Returns:
        Process exit status.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    if args.action == "status":
        path = args.output_dir / "status.json"
        print(json.dumps(json.loads(path.read_text(encoding="utf-8")), indent=2))
        return 0
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    overrides: dict[str, str] = {}
    for item in args.env:
        key, separator, value = item.partition("=")
        if not separator:
            raise ValueError(f"--env requires KEY=VALUE: {item!r}")
        overrides[key] = value
    script = launch_job(
        args.session,
        args.output_dir,
        args.cwd,
        command,
        env_overrides=overrides,
        input_paths=args.input,
    )
    LOGGER.info("launched session %s; run script: %s", args.session, script)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
