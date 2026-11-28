"""Readable plots for a frozen-data estimator discretization diagnostic."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_comparison(
    s: np.ndarray, curves: dict, center_curves: dict, output: Path
) -> None:
    """Write unsmoothed full-path/early plots and separate cell-center diagnostic.

    Args:
        s: Frozen projection-bin centers.
        curves: MBAR and WHAM curves projected at original frame values.
        center_curves: Separately labeled occupied-cell-center projections.
        output: Existing result directory receiving PNG and PDF figures.
    """
    plt.rcParams.update({"font.size": 16, "axes.labelsize": 18, "axes.titlesize": 19})
    colors = ["black", "#cf5543", "#d69d00", "#39995b", "#3578b5", "#9257ad"]
    for kind, selected in [
        ("original_frames", curves),
        (
            "cell_centers",
            {next(iter(curves)): next(iter(curves.values())), **center_curves},
        ),
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(17, 7), constrained_layout=True)
        for ax in axes:
            for i, (name, values) in enumerate(selected.items()):
                ax.plot(
                    s,
                    values,
                    color=colors[i],
                    lw=2.4 if i == 0 else 1.8,
                    ls="--" if "shift" in name else "-",
                    marker="o",
                    ms=3.5,
                    label=name.replace(" A", " Å"),
                )
            ax.set_xlabel("Hard path projection s")
            ax.set_ylabel("Relative PMF (kcal mol$^{-1}$)")
            ax.grid(alpha=0.2)
        axes[0].set_title("Whole occupied path")
        axes[1].set_title("Initial region: 0 ≤ s ≤ 0.3")
        axes[1].set_xlim(0, 0.3)
        early = (s >= 0) & (s <= 0.3)
        finite = np.concatenate(
            [v[early][np.isfinite(v[early])] for v in selected.values()]
        )
        axes[1].set_ylim(finite.min() - 0.8, finite.max() + 0.8)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=15)
        fig.suptitle(
            "Finite-bin WHAM vs MBAR: "
            + (
                "original-frame projection"
                if kind == "original_frames"
                else "separate cell-center projection diagnostic"
            ),
            fontsize=20,
        )
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"comparison_{kind}.{suffix}", dpi=180)
        plt.close(fig)
