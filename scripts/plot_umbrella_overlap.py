"""Plot normalized umbrella histograms and adjacent-window overlap."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def main() -> None:
    """Render a two-panel umbrella-overlap diagnostic."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()
    samples = json.loads(args.samples.read_text())
    summary = json.loads(args.summary.read_text())
    edges = np.asarray(summary["bin_edges_A"])
    centers = 0.5 * (edges[:-1] + edges[1:])
    figure, (histogram_axis, overlap_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        gridspec_kw={"height_ratios": [3, 1]},
        layout="constrained",
    )
    colors = plt.cm.viridis(np.linspace(0, 1, len(samples)))
    for color, row in zip(colors, samples, strict=True):
        counts, _ = np.histogram(row["samples_A"], bins=edges, density=True)
        histogram_axis.step(centers, counts, where="mid", color=color, alpha=0.8)
    histogram_axis.set(
        ylabel="Probability density",
        title=args.title,
    )
    histogram_axis.grid(alpha=0.2)
    overlaps = [row["neighbor_overlap"] for row in summary["windows"][:-1]]
    labels = [row["window"] for row in summary["windows"][:-1]]
    overlap_axis.bar(range(len(overlaps)), overlaps, color=colors[:-1])
    overlap_axis.set_xticks(range(len(labels)), labels, rotation=60, ha="right")
    overlap_axis.set(
        xlabel="Window paired with its next neighbor",
        ylabel="Histogram\nintersection",
        ylim=(0, 1),
    )
    overlap_axis.grid(axis="y", alpha=0.2)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
