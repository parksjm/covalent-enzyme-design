#!/usr/bin/env python3
"""Prepare a bounded 30-to-60 ps continuation without launching native workers."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.relay_umbrella_extension import prepare_extension


def main() -> None:
    """Write reviewed continuation inputs and immutable runtime/source identities."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-plan", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--derived-repair-snapshot", required=True, type=Path)
    args = parser.parse_args()
    summary = prepare_extension(
        args.original_plan,
        args.output,
        Path(__file__).with_name("guard_relay_umbrella_extension.py"),
        args.derived_repair_snapshot,
    )
    logging.info(
        "Prepared %s windows, target60ps; no native launch", summary["windows"]
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
