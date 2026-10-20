"""Export completed Nyl12 AF3 TI1 models with exact atom names and audits."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.af3_postprocess import canonicalize_af3_ti1  # noqa: E402
from nyl12_covalent.chemistry import build_poly_xama_maxa, load_pa6  # noqa: E402

SOURCE = ROOT / (
    "output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/1_relaxation/"
    "pa6_l4z_l2_nyl12_wt_model_41"
)
OUTPUT_ROOT = ROOT / "output_nyl12_covalent_qmmm_20261019"
LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Canonicalize one saved AF3 TI1 prediction and record its geometry."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("pa6_l2", "poly_xama_maxa"))
    parser.add_argument("--cif", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    substrate = (
        load_pa6(SOURCE / "lig_sub/lig_h.mol2")
        if args.target == "pa6_l2" else build_poly_xama_maxa()
    )
    audit = canonicalize_af3_ti1(args.cif, substrate, args.output)
    LOGGER.info("Canonicalized %s: %s", args.target, audit["distances_A"])


if __name__ == "__main__":
    main()
