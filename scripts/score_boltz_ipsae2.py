#!/usr/bin/env python3
"""Plan, poll or watch separately configured native ipSAE2 sidecar cases."""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

from nyl12_covalent.ipsae2_sidecar import poll_cases, read_json


def main() -> None:
    """Run the optional scorer with an explicit plan and qualification gate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true")
    mode.add_argument("--once", action="store_true")
    mode.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=float, default=60.0)
    args = parser.parse_args()
    if args.interval <= 0:
        parser.error("--interval must be positive")
    logging.basicConfig(level=logging.INFO)
    config = read_json(args.config)
    while True:
        result = poll_cases(config, plan=args.plan)
        logging.info("%s", json.dumps(result, sort_keys=True))
        terminal = result.get("cases", {})
        all_terminal = bool(terminal) and all(
            case["status"] in {"completed", "failed"} for case in terminal.values()
        )
        if not args.watch or all_terminal:
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
