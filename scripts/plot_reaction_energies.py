"""Render the energy data of a verified reaction-path inspection snapshot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    """Plot relative total and driving-bias-removed energies in kcal/mol."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = json.loads(args.data.read_text())
    x = list(range(len(rows)))
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    ax.plot(
        x,
        [r["relative_total_kcal_mol"] for r in rows],
        "o-",
        linewidth=2,
        label="Including driving restraints",
    )
    ax.plot(
        x,
        [r["relative_minus_driving_kcal_mol"] for r in rows],
        "s--",
        label="Driving bias removed",
    )
    for i, row in enumerate(rows):
        if not row["accepted"]:
            ax.axvspan(i - 0.15, i + 0.15, color="red", alpha=0.15)
        if i and row["stage"] != rows[i - 1]["stage"]:
            ax.axvline(i - 0.5, color="black", linestyle=":", alpha=0.5)
    ax.set_xticks(x, [r["window"] for r in rows], rotation=55, ha="right")
    ax.set(
        xlabel="Reaction-coordinate window",
        ylabel="Relative energy (kcal/mol)",
        title="Asp272 relay: finalized endpoints (red shading = rejected window)",
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.suptitle(
        "First-window reference; MM positional restraints retained.\n"
        "Construction scan, not a validated barrier.",
        fontsize=10,
    )
    fig.savefig(args.output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
