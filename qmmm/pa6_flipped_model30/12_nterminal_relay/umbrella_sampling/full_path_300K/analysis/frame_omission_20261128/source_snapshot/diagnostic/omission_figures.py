"""Curiosity-only PMF sensitivity figures with fixed baseline reference."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_profiles(s: np.ndarray, curves: dict, output: Path) -> None:
    """Plot whole and initial profiles of the three precisely labeled estimators.

    Args:
        s: Frozen projected-coordinate bin centers.
        curves: Baseline, full-refit and fixed-original-weight profile dictionaries.
        output: New diagnostic result directory.
    """
    plt.rcParams.update({"font.size": 17, "axes.labelsize": 19, "legend.fontsize": 16})
    fig, axes = plt.subplots(1, 2, figsize=(17, 7), constrained_layout=True)
    settings = [
        ("baseline", "Original MBAR · 82,000 frames", "#222222", "-", 2.6),
        ("full_refit", "One frame omitted · full MBAR refit", "#c62828", "-", 2.4),
        ("fixed_weights", "Omitted · fixed original weights", "#1976a3", "--", 2.1),
    ]
    for ax in axes:
        for key, label, color, style, width in settings:
            ax.plot(s, curves[key]["pmf"], label=label, color=color, ls=style, lw=width)
        ax.set(
            xlabel="Projected coordinate, s", ylabel="Relative PMF (kcal mol$^{-1}$)"
        )
        ax.grid(alpha=0.18)
    axes[0].set_title("Whole projected profile")
    axes[1].set_title("Initial region · affected bin at s = 0.2625")
    axes[1].set_xlim(-0.025, 0.325)
    use = (s >= -0.025) & (s <= 0.325)
    values = np.concatenate([v["pmf"][use] for v in curves.values()])
    axes[1].set_ylim(np.nanmin(values) - 0.7, np.nanmax(values) + 0.7)
    axes[1].axvline(0.2625, color="#c62828", alpha=0.35, lw=1)
    axes[0].legend(loc="lower right", frameon=False)
    fig.suptitle(
        "DIAGNOSTIC ONLY · NOT FOR PUBLICATION\n"
        "Single-observation omission; all curves zero at s = 0.0125",
        fontsize=19,
    )
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"single_frame_omission.{suffix}", dpi=180)
    plt.close(fig)
