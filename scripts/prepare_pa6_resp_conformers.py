"""Prepare explicit PA6 L4 zwitterion conformers for RESP geometry sampling."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_resp_preparation import (  # noqa: E402
    prepare_pa6_resp_conformers,
)


def _conformer(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("conformer must be LABEL=/absolute/path.pdb")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("conformer must contain a label and path")
    path = Path(raw_path)
    if not path.is_absolute():
        raise argparse.ArgumentTypeError("conformer PDB path must be absolute")
    return label, path


def main() -> None:
    """Parse explicit chemistry, conformers, constraints, and output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-mol2", type=Path, required=True)
    parser.add_argument(
        "--conformer",
        action="append",
        type=_conformer,
        required=True,
        help="Repeat LABEL=/absolute/path.pdb for every declared conformer",
    )
    parser.add_argument("--torsions-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    conformers = dict(args.conformer)
    if len(conformers) != len(args.conformer):
        parser.error("conformer labels must be unique")
    logging.basicConfig(level=logging.INFO)
    prepare_pa6_resp_conformers(
        args.source_mol2,
        conformers,
        args.torsions_json,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
