"""Create a provenance-checked MCCE protein input staging directory."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.mcce_preparation import prepare_input  # noqa: E402

LOGGER = logging.getLogger(__name__)


def main(argv: Sequence[str] | None = None) -> None:
    """Validate CLI arguments and prepare a new MCCE input directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument(
        "--termini",
        required=True,
        type=Path,
        help='JSON mapping chains to residue numbers, e.g. {"A":{"n":1,"c":232}}',
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    termini = json.loads(args.termini.read_text(encoding="utf-8"))
    report = prepare_input(args.source, args.output, args.source_sha256, termini)
    LOGGER.info(
        "Prepared %s (%d atom names normalized)",
        args.output,
        len(report["changed_atom_names"]),
    )


if __name__ == "__main__":
    main()
