#!/usr/bin/env python3
"""Calculate pointwise statistical errors for a frozen relay PMF checkpoint."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.relay_pmf_uncertainty import report_uncertainty
from nyl12_covalent.relay_umbrella_analysis_io import digest, load_block


def main() -> None:
    """Validate immutable inputs and write a separate, source-bound report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=0)
    parser.add_argument("--bootstrap-workers", type=int, default=1)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve existing analysis; choose a new output")
    config = json.loads(args.request.read_text())
    report = json.loads(args.report.read_text())
    if digest(args.request) != report["request_sha256"]:
        raise ValueError("Frozen request hash does not match source report")
    blocks = [
        load_block(
            (args.request.parent / receipt).resolve(strict=True),
            temperature=config["temperature_K"],
            hamiltonian=config["physical_hamiltonian_sha256"],
            diagnostic_columns=config["diagnostic_columns"],
        )
        for receipt in config["accepted_block_receipts"]
    ]
    result = report_uncertainty(
        config,
        blocks,
        report,
        bootstrap_replicates=args.bootstrap_replicates,
        bootstrap_workers=args.bootstrap_workers,
    )
    result.update(
        report_sha256=digest(args.report),
        request_sha256=digest(args.request),
        report_source=str(args.report.resolve()),
        source_sha256={
            str(path): digest(path)
            for path in (
                Path(__file__),
                Path(__file__).resolve().parents[1]
                / "src/nyl12_covalent/relay_pmf_uncertainty.py",
            )
        },
    )
    args.output.parent.mkdir(exist_ok=True, parents=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    logging.info("Wrote conditional statistical errors: %s", args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
