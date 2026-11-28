#!/usr/bin/env python3
"""Prepare six completed native cases; retain fresh mappings and pairing evidence."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.ipsae2_prepare import prepare_completed_campaign


def main() -> None:
    """Prepare explicit production sidecar inputs without launching scoring."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "campaign",
        "output",
        "runtime-path",
        "qualification-path",
        "pairing-adapter",
        "boltz-python",
        "boltz-source-root",
        "root",
    ):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    config = prepare_completed_campaign(
        **{key: value.resolve() for key, value in vars(args).items()}
    )
    logging.info("Prepared sidecar configuration: %s", config)


if __name__ == "__main__":
    main()
