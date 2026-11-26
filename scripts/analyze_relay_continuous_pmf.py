#!/usr/bin/env python3
"""Reanalyze a frozen PMF checkpoint with one production burn-in per window."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.relay_continuous_pmf import analyze_continuous
from nyl12_covalent.relay_umbrella_analysis_io import digest, load_block


def main() -> None:
    """Verify sources and publish separate analysis/uncertainty JSON files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--parent-report", required=True, type=Path)
    parser.add_argument("--output-stem", required=True, type=Path)
    parser.add_argument("--replicates", required=True, type=int)
    parser.add_argument("--workers", required=True, type=int)
    args = parser.parse_args()
    report_path = args.output_stem.with_name(args.output_stem.name + "_analysis.json")
    uncertainty_path = args.output_stem.with_name(
        args.output_stem.name + "_uncertainty.json"
    )
    if report_path.exists() or uncertainty_path.exists():
        raise FileExistsError("Retain previous analyses; choose a new output stem")
    parent = json.loads(args.parent_report.read_text())
    if digest(args.request) != parent["request_sha256"]:
        raise ValueError("Request does not match frozen parent report")
    config = json.loads(args.request.read_text())
    blocks = [
        load_block(
            (args.request.parent / filename).resolve(strict=True),
            temperature=config["temperature_K"],
            hamiltonian=config["physical_hamiltonian_sha256"],
            diagnostic_columns=config["diagnostic_columns"],
        )
        for filename in config["accepted_block_receipts"]
    ]
    report, uncertainty = analyze_continuous(
        config,
        blocks,
        parent,
        replicates=args.replicates,
        workers=args.workers,
    )
    metadata = dict(
        parent_report_sha256=digest(args.parent_report),
        request_sha256=digest(args.request),
        source_sha256={
            str(path): digest(path)
            for path in [
                Path(__file__),
                Path(__file__).resolve().parents[1]
                / "src/nyl12_covalent/relay_continuous_pmf.py",
            ]
        },
    )
    report.update(metadata)
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    uncertainty.update(metadata, report_sha256=digest(report_path))
    uncertainty_path.write_text(
        json.dumps(uncertainty, indent=2, allow_nan=False) + "\n"
    )
    logging.info("Continuous-window PMF: %s", report_path)
    logging.info("Selected configurations: %s", sum(report["N_k"]))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
