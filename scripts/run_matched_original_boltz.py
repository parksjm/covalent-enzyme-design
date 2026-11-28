"""Run explicit matched precursor search, staging, audit, or stock prediction."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.matched_original_boltz import journal, prepare, search


def main() -> None:
    """Execute only the selected stage for one explicitly configured candidate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--stage", choices=["search", "prepare", "predict"], required=True
    )
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()
    if args.seed is not None and args.stage != "predict":
        parser.error("--seed is only valid for the predict stage")
    plan = json.loads(args.config.read_text())
    output = args.output.resolve()
    journal("intent", f"Matched Boltz stage {args.stage}: {args.candidate}")
    if args.stage == "search":
        search(plan, output, args.candidate)
    elif args.stage == "prepare":
        prepare(plan, output, args.candidate)
    else:
        from nyl12_covalent.matched_original_native import predict_candidate

        predict_candidate(plan, output, args.candidate, seed=args.seed)
    journal("outcome", f"Matched Boltz stage {args.stage} completed: {args.candidate}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
