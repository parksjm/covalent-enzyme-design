"""Generate a new PA6 constrained Gaussian geometry attempt without running it."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.pa6_resp_restart import (  # noqa: E402
    GaussianSettings,
    prepare_pa6_resp_restart,
)


def main() -> None:
    """Parse all required method and resource settings and write new inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--basis", required=True)
    parser.add_argument("--solvent", required=True, help="Explicitly as SMD:Water")
    parser.add_argument("--max-cycles", type=int, required=True)
    parser.add_argument("--cores", type=int, required=True)
    parser.add_argument("--memory-gb", type=int, required=True)
    parser.add_argument(
        "--calc-fc", action=argparse.BooleanOptionalAction, required=True
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    settings = GaussianSettings(
        method=args.method,
        basis=args.basis,
        solvent=args.solvent,
        max_cycles=args.max_cycles,
        cores=args.cores,
        memory_gb=args.memory_gb,
        calc_fc=args.calc_fc,
    )
    prepare_pa6_resp_restart(args.source_manifest, args.output_dir, settings)


if __name__ == "__main__":
    main()
