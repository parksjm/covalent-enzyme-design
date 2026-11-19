#!/usr/bin/env python3
"""Build the source-pinned solvated Amber system for PA6 Model 23."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.es_system_preparation import prepare_es_system  # noqa: E402


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdb", type=Path, required=True)
    parser.add_argument("--ligand-library", type=Path, required=True)
    parser.add_argument("--ligand-frcmod", type=Path, required=True)
    parser.add_argument("--tleap", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-label", default="Model 23")
    parser.add_argument("--source-snapshot-name", default="source_model_23.pdb")
    return parser.parse_args()


def main() -> int:
    """Prepare the system and log its audited summary."""
    args = parse_args()
    report = prepare_es_system(
        args.source_pdb,
        args.ligand_library,
        args.ligand_frcmod,
        args.tleap,
        args.output,
        args.source_label,
        args.source_snapshot_name,
    )
    LOGGER.info("Prepared Nyl12 E-S system: %s", json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
