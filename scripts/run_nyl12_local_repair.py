"""Run a separately recorded alpha/beta TI1 restrained local repair in tmux."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.local_repair_run import run_repair  # noqa: E402


def main() -> None:
    """Run the explicitly requested, fresh local-repair calculation."""
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("topology", "coordinates", "cif", "amber-bin", "output"):
        parser.add_argument(f"--{argument}", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_repair(
        args.topology,
        args.coordinates,
        args.cif,
        args.amber_bin,
        args.output.absolute(),
        args.prepare_only,
    )


if __name__ == "__main__":
    main()
