#!/usr/bin/env python3
"""Convert one historical Nyl12 whole-PA6 complex to modular residues."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.modular_pa6_pdb import convert_legacy_pa6_l4_pdb


def main() -> int:
    """Parse inputs and write an immutable modular-complex conversion."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdb", required=True, type=Path)
    parser.add_argument("--source-mol2", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--release-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    report = convert_legacy_pa6_l4_pdb(
        args.source_pdb,
        args.source_mol2,
        args.mapping,
        args.release_dir,
        args.output_dir,
    )
    logging.info("Converted %d PA6 atoms", len(report["atom_mapping"]))
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
