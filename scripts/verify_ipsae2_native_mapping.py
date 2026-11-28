#!/usr/bin/env python3
"""Write a strict native mapping receipt using the existing Boltz environment."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.ipsae2_native_mapping import qualify_native_mapping


def main() -> None:
    """Validate explicit native artifacts and persist a new receipt."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "cif",
        "pae",
        "confidence",
        "plddt",
        "structure",
        "boltz-source-root",
        "output",
    ):
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args()
    receipt = qualify_native_mapping(
        cif_path=args.cif,
        pae_path=args.pae,
        confidence_path=args.confidence,
        plddt_path=args.plddt,
        structure_path=args.structure,
        boltz_source_root=args.boltz_source_root,
    )
    with args.output.open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    logging.basicConfig(level=logging.INFO)
    logging.info(
        "Native mapping accepted=%s; receipt=%s", receipt["accepted"], args.output
    )
    if not receipt["accepted"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
