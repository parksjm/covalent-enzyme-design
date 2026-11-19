#!/usr/bin/env python3
"""Align a canonical static Nyl12 E-S PDB to the catalytic reference."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.static_es_alignment import align_static_es_model  # noqa: E402


def main() -> None:
    """Run the static E-S catalytic alignment from explicit paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_cif", type=Path)
    parser.add_argument("model_pdb", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--model-chain", default="A")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)
    align_static_es_model(
        args.reference_cif,
        args.model_pdb,
        args.output_dir,
        args.model_chain,
    )


if __name__ == "__main__":
    main()
