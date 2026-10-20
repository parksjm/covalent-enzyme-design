"""Audit the Nyl12 TI1 QM atom selection without launching an engine."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.qm_region import audit_qm_region  # noqa: E402


def main() -> None:
    """Parse audit inputs and write the structural selection report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("topology", type=Path, help="Amber prmtop to inspect")
    parser.add_argument(
        "--source",
        type=Path,
        action="append",
        default=[],
        help="Additional source file to hash; may be repeated",
    )
    parser.add_argument("--output", type=Path, help="Write JSON report to this path")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    report = audit_qm_region(args.topology, args.source)
    serialized = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(serialized, encoding="utf-8")
    else:
        sys.stdout.write(serialized)


if __name__ == "__main__":
    main()
