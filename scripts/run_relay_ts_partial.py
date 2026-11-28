"""Run the bounded partial-QM relay TS curvature diagnostic."""

import argparse
import logging
from pathlib import Path

from nyl12_covalent.relay_ts_partial import partial_screen


def main() -> None:
    """Parse explicit immutable inputs and execute the diagnostic."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    partial_screen(
        args.model.resolve(),
        args.restart.resolve(),
        args.output.resolve(),
        args.workers,
    )


if __name__ == "__main__":
    main()
