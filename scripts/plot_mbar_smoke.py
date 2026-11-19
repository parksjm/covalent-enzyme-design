#!/usr/bin/env python3
"""Plot a clearly labelled provisional PMF from an MBAR smoke test."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    """Render provisional PMF values and analytical uncertainties."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--title", required=True)
    parser.add_argument(
        "--banner",
        default="PIPELINE SMOKE TEST — NOT CONVERGED",
        help="Explicit interpretation banner shown on the plot",
    )
    args = parser.parse_args()
    result = json.loads(args.summary.read_text())
    coordinate = [row["coordinate_A"] for row in result["provisional_pmf"]]
    energy = [row["free_energy_kcal_mol"] for row in result["provisional_pmf"]]
    uncertainty = [row["uncertainty_kcal_mol"] for row in result["provisional_pmf"]]
    figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    axis.errorbar(coordinate, energy, yerr=uncertainty, marker="o", ms=3, lw=1)
    axis.set(
        xlabel="Umbrella coordinate (Å)",
        ylabel="Provisional free energy (kcal/mol)",
        title=args.title,
    )
    axis.text(
        0.02,
        0.98,
        args.banner,
        transform=axis.transAxes,
        va="top",
        color="firebrick",
        weight="bold",
    )
    axis.grid(alpha=0.25)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
