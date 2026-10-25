"""Render declared reaction distances using the explicit plotting environment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    """Plot endpoint distances and mark rejected windows and jump flags."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    args = parser.parse_args()
    data = json.loads(args.data.read_text())
    rows = data["rows"]
    groups = list(dict.fromkeys(pair["group"] for pair in data["pairs"].values()))
    fig, axes = plt.subplots(
        len(groups),
        1,
        figsize=(10, 2.8 * len(groups)),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )
    x = list(range(len(rows)))
    for group, ax in zip(groups, axes[:, 0], strict=True):
        for name, definition in data["pairs"].items():
            if definition["group"] == group:
                ax.plot(
                    x,
                    [r["distances_A"][name] for r in rows],
                    "o-",
                    markersize=4,
                    label=name,
                )
        for index, row in enumerate(rows):
            if not row["accepted"]:
                ax.axvspan(index - 0.15, index + 0.15, color="red", alpha=0.12)
            if index and row["stage"] != rows[index - 1]["stage"]:
                ax.axvline(index - 0.5, color="black", linestyle="--", alpha=0.5)
        flagged = {
            f["to"]
            for f in data["jumps"]
            if data["pairs"][f["distance"]]["group"] == group
        }
        for index, row in enumerate(rows):
            if row["window"] in flagged:
                ax.axvline(index, color="darkorange", linestyle=":", alpha=0.8)
        ax.set(ylabel="Distance (Å)", title=group)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8, ncol=3)
    axes[-1, 0].set_xticks(x, [r["window"] for r in rows], rotation=55, ha="right")
    axes[-1, 0].set_xlabel(
        "Reaction-coordinate window (one optimized endpoint per step)"
    )
    fig.suptitle(
        f"Reaction distances | {len(data['jumps'])} changes > "
        f"{data['jump_threshold_A']:.2f} Å flagged\n"
        "Red shading: rejected endpoint; orange dotted: flagged distance change",
        fontsize=11,
    )
    fig.savefig(args.data.with_suffix(".png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
