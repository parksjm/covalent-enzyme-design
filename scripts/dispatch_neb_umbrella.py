#!/usr/bin/env python3
"""Prepare immutable two-host dispatch, or execute an explicit host assignment."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_umbrella_dispatch import launch, prepare, worker  # noqa: E402


def main() -> None:
    """Require explicit plan preparation and host-specific execution actions."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "launch", "worker"))
    parser.add_argument("--config", type=Path)
    parser.add_argument("--hosts", type=Path)
    parser.add_argument("--local-host", default="violet")
    parser.add_argument("--remote-host", default="azalea")
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--host")
    parser.add_argument("--worker", type=int)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.action == "prepare":
        if args.config is None or args.hosts is None:
            parser.error("prepare requires --config and --hosts")
        path = prepare(
            args.config.resolve(strict=True),
            json.loads(args.hosts.read_text()),
            args.local_host,
            args.remote_host,
        )
        logging.info("Prepared immutable dispatch plan: %s", path)
    else:
        if args.plan is None or args.host is None:
            parser.error("execution requires --plan and --host")
        path = args.plan.resolve(strict=True)
        if args.action == "launch":
            launch(path, args.host)
        elif args.worker is None:
            parser.error("worker requires --worker")
        else:
            worker(path, args.host, args.worker)


if __name__ == "__main__":
    main()
