#!/usr/bin/env python3
"""Run a provenance-recorded multi-block PyMBAR umbrella analysis."""

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
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.umbrella_mbar import (  # noqa: E402
    load_umbrella_sample_blocks,
    run_mbar_smoke,
)

LOGGER = logging.getLogger(__name__)


def _profile_change(
    previous: list[dict[str, float]], current: list[dict[str, float]]
) -> dict[str, float]:
    """Compare two zero-referenced PMFs on their common coordinate domain."""
    previous_x = np.asarray([row["coordinate_A"] for row in previous])
    previous_y = np.asarray([row["free_energy_kcal_mol"] for row in previous])
    current_x = np.asarray([row["coordinate_A"] for row in current])
    current_y = np.asarray([row["free_energy_kcal_mol"] for row in current])
    lower = max(float(previous_x.min()), float(current_x.min()))
    upper = min(float(previous_x.max()), float(current_x.max()))
    grid = current_x[(current_x >= lower) & (current_x <= upper)]
    if grid.size < 4:
        raise ValueError("PMFs have fewer than four common comparison points")
    old = np.interp(grid, previous_x, previous_y)
    new = np.interp(grid, current_x, current_y)
    old -= old.min()
    new -= new.min()
    difference = new - old
    return {
        "common_point_count": int(grid.size),
        "rms_change_kcal_mol": float(np.sqrt(np.mean(np.square(difference)))),
        "mean_absolute_change_kcal_mol": float(np.mean(np.abs(difference))),
        "maximum_absolute_change_kcal_mol": float(np.max(np.abs(difference))),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a nonempty list of uniform dictionaries as CSV."""
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(
    sources: list[tuple[Path, str]],
    output: Path,
    bins: int,
    plot_python: Path,
    title: str,
) -> None:
    """Analyze cumulative complete umbrella blocks and write diagnostics.

    Args:
        sources: Ordered campaign and exact-stage pairs.
        output: New analysis result directory.
        bins: Histogram bin count for the FES estimator.
        plot_python: Python executable containing Matplotlib.
        title: Human-readable plot title.
    """
    output.mkdir(parents=True, exist_ok=False)
    cumulative: list[dict[str, Any]] = []
    previous_profile: list[dict[str, float]] | None = None
    final: dict[str, Any] | None = None
    for source_count in range(1, len(sources) + 1):
        inputs, windows = load_umbrella_sample_blocks(
            sources[:source_count], require_successful_process=True
        )
        result = run_mbar_smoke(
            windows,
            temperature_k=float(inputs["temperature_K"]),
            force_constant_kcal_mol_a2=float(inputs["force_constant_kcal_mol_A2"]),
            bins=bins,
        )
        result["interpretation"] = (
            "interim multi-block PMF; convergence requires additional independent "
            "sampling and stability across cumulative blocks"
        )
        entry: dict[str, Any] = {
            "source_count": source_count,
            "raw_sample_count": result["raw_sample_count"],
            "decorrelated_sample_count": result["decorrelated_sample_count"],
            "mbar_overlap_scalar": result["mbar_overlap_scalar"],
            "minimum_bidirectional_adjacent_overlap": result[
                "minimum_bidirectional_adjacent_overlap"
            ],
            "pmf_range_kcal_mol": max(
                row["free_energy_kcal_mol"] for row in result["provisional_pmf"]
            ),
        }
        if previous_profile is not None:
            entry["change_from_previous"] = _profile_change(
                previous_profile, result["provisional_pmf"]
            )
        cumulative.append(entry)
        previous_profile = result["provisional_pmf"]
        final = result
    if final is None:
        raise ValueError("At least one source block is required")

    final["sources"] = [
        {
            "campaign": str(campaign.resolve()),
            "stage": stage,
            "inputs_sha256": hashlib.sha256(
                (campaign / "inputs.json").read_bytes()
            ).hexdigest(),
        }
        for campaign, stage in sources
    ]
    final["cumulative_convergence"] = cumulative
    summary = output / "mbar_summary.json"
    summary.write_text(json.dumps(final, indent=2) + "\n", encoding="utf-8")
    profile = output / "pmf.csv"
    _write_csv(profile, final["provisional_pmf"])
    convergence = output / "cumulative_convergence.json"
    convergence.write_text(json.dumps(cumulative, indent=2) + "\n", encoding="utf-8")
    overlap = output / "mbar_overlap_matrix.csv"
    labels = [row["window"] for row in final["windows"]]
    with overlap.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["source_window", *labels])
        for label, row in zip(labels, final["mbar_overlap_matrix"], strict=True):
            writer.writerow([label, *row])

    plot = output / "pmf.png"
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    subprocess.run(
        [
            str(plot_python),
            str(ROOT / "scripts/plot_mbar_smoke.py"),
            str(summary),
            str(plot),
            "--title",
            title,
            "--banner",
            "INTERIM ANALYSIS — NOT CONVERGED",
        ],
        check=True,
        env=environment,
    )
    artifacts = (summary, profile, convergence, overlap, plot)
    (output / "MANIFEST.sha256").write_text(
        "".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n"
            for path in artifacts
        ),
        encoding="utf-8",
    )
    LOGGER.info("Wrote multi-block MBAR analysis to %s", output)


def main() -> None:
    """Parse source blocks and run PyMBAR analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        nargs=2,
        action="append",
        metavar=("CAMPAIGN", "STAGE"),
        required=True,
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--plot-python", required=True, type=Path)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    sources = [(Path(campaign), stage) for campaign, stage in args.source]
    analyze(sources, args.output, args.bins, args.plot_python, args.title)


if __name__ == "__main__":
    main()
