"""Prepare or externally gate the exact eight-model unguided predictor comparison."""

import argparse
import logging
from pathlib import Path

from nyl12_covalent.relay_broad_predictions import (
    DEFAULT_ATOMS,
    DEFAULT_PROTECTION,
    prepare_comparison,
    run_comparison,
)


def main() -> None:
    """Prepare inputs by default; explicit --run is reserved for the gated controller."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--design-json", type=Path, required=True)
    parser.add_argument("--protection", type=Path, default=DEFAULT_PROTECTION)
    parser.add_argument("--native-atom-table", type=Path, default=DEFAULT_ATOMS)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    manifest = prepare_comparison(
        args.campaign, args.design_json, args.protection, args.native_atom_table
    )
    if args.run:
        run_comparison(manifest, args.campaign)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
