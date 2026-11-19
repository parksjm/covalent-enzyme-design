#!/usr/bin/env python3
"""Transfer an accepted PA6 DFT geometry to its RESP/GAFF2 MOL2."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_validation_geometry import (  # noqa: E402
    prepare_pa6_validation_geometry,
)

LOGGER = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Run accepted-coordinate transfer.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-mol2", required=True, type=Path)
    parser.add_argument("--typed-mol2", required=True, type=Path)
    parser.add_argument("--gaussian-log", required=True, type=Path)
    parser.add_argument("--acceptance-json", required=True, type=Path)
    parser.add_argument("--output-mol2", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    prepare_pa6_validation_geometry(
        args.source_mol2,
        args.typed_mol2,
        args.gaussian_log,
        args.acceptance_json,
        args.output_mol2,
        args.manifest,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
