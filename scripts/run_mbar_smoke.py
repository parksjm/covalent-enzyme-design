#!/usr/bin/env python3
"""Run a provenance-recorded, explicitly provisional MBAR pipeline test."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.umbrella_mbar import (  # noqa: E402
    load_umbrella_samples,
    run_mbar_smoke,
)

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Parse arguments, run MBAR diagnostics, and write hashed outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--active-snapshot", action="store_true")
    parser.add_argument("--plot-python", type=Path, required=True)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    inputs, windows = load_umbrella_samples(
        args.campaign,
        args.stage,
        require_successful_process=not args.active_snapshot,
    )
    result = run_mbar_smoke(
        windows,
        temperature_k=float(inputs["temperature_K"]),
        force_constant_kcal_mol_a2=float(inputs["force_constant_kcal_mol_A2"]),
        bins=args.bins,
    )
    result["campaign"] = str(args.campaign.resolve())
    result["stage"] = args.stage
    result["active_snapshot"] = args.active_snapshot
    result["inputs_sha256"] = hashlib.sha256(
        (args.campaign / "inputs.json").read_bytes()
    ).hexdigest()

    args.output.mkdir(parents=True, exist_ok=False)
    summary = args.output / "mbar_smoke_summary.json"
    summary.write_text(json.dumps(result, indent=2) + "\n")
    profile = args.output / "provisional_pmf.csv"
    with profile.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result["provisional_pmf"][0]))
        writer.writeheader()
        writer.writerows(result["provisional_pmf"])
    overlap = args.output / "mbar_overlap_matrix.csv"
    with overlap.open("w", newline="") as handle:
        writer = csv.writer(handle)
        labels = [row["window"] for row in result["windows"]]
        writer.writerow(["source_window", *labels])
        for label, row in zip(labels, result["mbar_overlap_matrix"], strict=True):
            writer.writerow([label, *row])

    plot = args.output / "provisional_pmf.png"
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    subprocess.run(
        [
            str(args.plot_python),
            str(Path(__file__).with_name("plot_mbar_smoke.py")),
            str(summary),
            str(plot),
            "--title",
            args.title,
        ],
        check=True,
        env=environment,
    )
    outputs = (summary, profile, overlap, plot)
    manifest = args.output / "MANIFEST.sha256"
    manifest.write_text(
        "".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n"
            for path in outputs
        )
    )
    LOGGER.info("Wrote provisional MBAR smoke test to %s", args.output)


if __name__ == "__main__":
    main()
