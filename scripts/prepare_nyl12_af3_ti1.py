"""Prepare product-gated AF3 jobs for Nyl12 tetrahedral covalent adducts."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.af3 import prepare_af3_ti1  # noqa: E402
from nyl12_covalent.chemistry import (  # noqa: E402
    build_poly_xama_maxa,
    load_pa6,
)

SOURCE = ROOT / (
    "output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/1_relaxation/"
    "pa6_l4z_l2_nyl12_wt_model_41"
)
MSA_SOURCE = Path(
    "/data/jpc/code/alphafold3/nyl12_polyamide_altlength_stoich_20261009/"
    "af3_output_pa6_l4z_single/"
    "nyl12_wt_pa6_l4z_single_unrestrained_af3_model_03/"
    "nyl12_wt_pa6_l4z_single_unrestrained_af3_model_03_data.json"
)


def main() -> None:
    """Write a custom CCD, covalent AF3 input, runner, and provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("pa6_l2", "poly_xama_maxa"))
    parser.add_argument("--pose", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--msa-source", type=Path, default=MSA_SOURCE)
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--gpu", type=int, default=0)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.target == "pa6_l2":
        substrate = load_pa6(SOURCE / "lig_sub/lig_h.mol2")
        pose = args.pose or SOURCE / "complex.pdb"
    else:
        substrate = build_poly_xama_maxa()
        if args.pose is None:
            parser.error("--pose is required for poly(X,AMA)")
        pose = args.pose
    output = args.output or (
        ROOT / "output_nyl12_covalent_qmmm_20261019" / f"{args.target}_af3_ti1"
    )
    prepare_af3_ti1(substrate, pose, args.msa_source, output, args.seed, args.gpu)
    logging.info("Wrote AF3 tetrahedral job: %s", output)


if __name__ == "__main__":
    main()
