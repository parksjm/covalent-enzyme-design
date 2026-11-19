#!/usr/bin/env python3
"""Prepare the template-conditioned Nyl12 PA6 E-S Boltz pilot inputs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from nyl12_covalent.template_boltz import prepare_pilot  # noqa: E402

DEFAULT_SOURCE = Path(
    "/data/jpc/enzyme-substrate-modeling/"
    "output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/0_boltz/input/"
    "yaml_incremental_msa/pa6_l4z_l2_nyl12_wt.yaml"
)
DEFAULT_TEMPLATE = Path("/data/jpc/covalent-enzyme-design/cif/n12c.cif")
DEFAULT_OUTPUT = Path(
    "/data/jpc/covalent-enzyme-design/inputs/template_conditioned_es_pilot"
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-yaml", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--template-cif", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed-base", type=int, default=2026102700)
    return parser.parse_args()


def main() -> None:
    """Prepare compact inputs without launching prediction jobs."""
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    prepare_pilot(
        source_yaml=args.source_yaml,
        template_cif=args.template_cif,
        output_dir=args.output_dir,
        seed_base=args.seed_base,
    )


if __name__ == "__main__":
    main()
