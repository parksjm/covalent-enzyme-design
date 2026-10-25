"""Prepare or explicitly run the force-gated, 25-cycle Nyl12 ES relaxation."""

from __future__ import annotations
import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.es_dftb_relaxation import (  # noqa: E402
    prepare_es_relaxation,
    run_prepared_relaxation,
)


def main() -> None:
    """Parse explicit force-gated sources; engine execution requires --run."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["engine-config", "force-check", "region-config", "sander", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument(
        "--run",
        action="store_true",
        help="Execute only after all gates pass; default prepares only.",
    )
    parser.add_argument(
        "--pilot-decision",
        type=Path,
        help="Explicit bounded analytic-pilot decision; default requires accepted FD proof.",
    )
    args = parser.parse_args()
    prepared = prepare_es_relaxation(
        args.engine_config,
        args.force_check,
        args.region_config,
        args.sander,
        args.output,
        args.pilot_decision,
    )
    if args.run:
        run_prepared_relaxation(prepared)
    logging.getLogger(__name__).info(
        "Prepared bounded ES pilot; engine run requested=%s", args.run
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
