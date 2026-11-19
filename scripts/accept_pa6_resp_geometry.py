"""Gate one PA6 DFT geometry and prepare its fixed-geometry HF/MK ESP input."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_resp_acceptance import prepare_pa6_hf_esp  # noqa: E402


def main() -> None:
    """Parse explicit source, conformer, threshold, and output arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-mol2", type=Path, required=True)
    parser.add_argument("--preparation-manifest", type=Path, required=True)
    parser.add_argument("--conformer", required=True)
    parser.add_argument("--gaussian-log", type=Path, required=True)
    parser.add_argument("--torsion-tolerance-degrees", type=float, required=True)
    parser.add_argument("--acceptance-json", type=Path, required=True)
    parser.add_argument("--esp-input", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    prepare_pa6_hf_esp(
        args.source_mol2,
        args.preparation_manifest,
        args.conformer,
        args.gaussian_log,
        args.torsion_tolerance_degrees,
        args.acceptance_json,
        args.esp_input,
    )


if __name__ == "__main__":
    main()
