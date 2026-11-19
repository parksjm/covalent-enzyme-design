#!/usr/bin/env python3
"""Fit shared PA6 L4 RESP charges and prepare non-certified GAFF2 candidates."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_resp_fit import fit_pa6_multiconformer_resp  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-mol2", required=True, type=Path)
    parser.add_argument("--esp", required=True, action="append", type=Path)
    parser.add_argument("--amber-bin", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the PA6 shared RESP fitting workflow.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Process exit status.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _parser().parse_args(argv)
    fit_pa6_multiconformer_resp(
        args.source_mol2, args.esp, args.amber_bin, args.output_dir
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, FileExistsError, RuntimeError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
