#!/usr/bin/env python3
"""Build and audit the reviewed modular-PA6 E-S QM region."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.modular_es_qm_region import build_modular_es_region  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Parse paths and write the audited modular QM-region configuration."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    result = build_modular_es_region(
        args.template.resolve(),
        args.topology.resolve(),
        args.coordinates.resolve(),
        args.output.resolve(),
    )
    LOGGER.info("Built modular E-S QM region: %s", result)


if __name__ == "__main__":
    main()
