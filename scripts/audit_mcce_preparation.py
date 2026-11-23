#!/usr/bin/env python3
"""Audit the scientific and process acceptance of an MCCE step1 run."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.mcce_validation import audit_preparation  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the step1 audit and write a new exclusive JSON report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--termini", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    expected_termini = json.loads(args.termini.read_text(encoding="utf-8"))
    report = audit_preparation(args.folder, args.input, expected_termini)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    LOGGER.info("Wrote acceptance report to %s", args.output)
    return 0 if report["accepted"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
