"""Plot five physical distances without connecting independent histories."""

from __future__ import annotations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from nyl12_covalent.reaction_distance_stats import DISTANCES, LABELS


def save(fig: plt.Figure, folder: Path, name: str) -> None:
    """Retain raster and vector versions of a diagnostic."""
    for extension in ("png", "pdf"):
        fig.savefig(folder / f"{name}.{extension}", dpi=170)
    plt.close(fig)


def plot_all(
    folder: Path,
    rows: list[dict],
    native: list[dict],
    data: dict,
    s: np.ndarray,
    centers: np.ndarray,
    statistics: list[dict],
) -> None:
    """Show raw samples, weighted quantiles, independent traces and local kinetics."""
    plt.rcParams.update(
        {"font.size": 12, "axes.spines.top": False, "axes.spines.right": False}
    )
    for zoom in (False, True):
        fig, axes = plt.subplots(
            5, 1, figsize=(13, 14), sharex=True, constrained_layout=True
        )
        for ax, key, label in zip(axes, DISTANCES, LABELS, strict=True):
            ax.scatter(
                s,
                data[key],
                s=2,
                alpha=0.045,
                color="0.3",
                rasterized=True,
                label="All retained frames (unweighted)",
            )
            values = [r for r in statistics if r["distance"] == key]
            x = np.array([r["s"] for r in values])
            ys = np.array(
                [
                    [r["weighted"][k] for k in ("p10_A", "median_A", "p90_A")]
                    for r in values
                ]
            )
            ax.fill_between(
                x,
                ys[:, 0],
                ys[:, 2],
                alpha=0.3,
                color="tab:orange",
                label="Weighted 10–90% distribution (not CI)",
            )
            ax.plot(x, ys[:, 1], color="tab:red", lw=1.8, label="Weighted median")
            ax.set(ylabel=f"{label}\n(Å)", ylim=(0, None))
            ax.grid(alpha=0.15)
        axes[0].legend(fontsize=10, ncol=1, loc="upper right")
        axes[-1].set(
            xlabel="Declared hard projected path progress, s",
            xlim=(0, 0.35) if zoom else (centers[0] - 0.0125, centers[-1] + 0.0125),
        )
        fig.suptitle(
            "Five direct reaction distances — frozen 100-history, 82,000-frame pool"
        )
        save(
            fig, folder, "five_distances_vs_s_early" if zoom else "five_distances_vs_s"
        )
    windows = list(
        dict.fromkeys(
            r["state_id"] for r in native if r["state_id"].startswith("main57/")
        )
    )
    fig, axes = plt.subplots(
        5, 1, figsize=(17, 13), sharex=True, constrained_layout=True
    )
    for window_index, window in enumerate(windows):
        group = [r for r in native if r["state_id"] == window]
        x = (
            window_index
            + np.array([r["production_time_ps"] for r in group]) / 70 * 0.90
        )
        for ax, key in zip(axes, DISTANCES, strict=True):
            ax.plot(
                x,
                [r[key] for r in group],
                lw=0.4,
                alpha=0.7,
                color=plt.get_cmap("tab10")(window_index % 10),
                rasterized=True,
            )
    for ax, label in zip(axes, LABELS, strict=True):
        ax.set_ylabel(f"{label}\n(Å)")
        ax.grid(alpha=0.12)
    axes[-1].set(
        xlabel="Main path window order; time increases separately within each window",
        xticks=np.arange(len(windows)) + 0.45,
        xticklabels=[
            w.split("/")[-1]
            .replace("activation_full_", "A")
            .replace("delivery_full_", "D")
            for w in windows
        ],
    )
    axes[-1].tick_params(axis="x", labelrotation=90, labelsize=6)
    fig.suptitle(
        "Main 57 windows, 0–70 ps each — independent traces; initial 5 ps included"
    )
    save(fig, folder, "five_distances_main57_70ps_independent_traces")
    plot_local(folder, native)


def plot_local(folder: Path, native: list[dict]) -> None:
    """Plot retained local histories with a readable delivery-axis label."""
    fig, axes = plt.subplots(
        6, 1, figsize=(13, 13), sharex=True, constrained_layout=True
    )
    for suffix, color in (
        ("3_02", "tab:blue"),
        ("3_03", "tab:red"),
        ("3_04", "tab:green"),
    ):
        group = [
            r
            for r in native
            if r["state_id"].endswith("activation_support_" + suffix)
            and 5 <= r["production_time_ps"] <= 10
        ]
        for ax, key in zip(axes, (*DISTANCES, "q_del_A"), strict=True):
            ax.plot(
                [r["production_time_ps"] for r in group],
                [r[key] for r in group],
                color=color,
                label="activation_support_" + suffix,
                lw=1,
            )
            ax.axvline(7.6, color="0.2", ls=":", lw=1)
    for ax, label in zip(axes, (*LABELS, "Delivery coordinate, q_del"), strict=True):
        ax.set_ylabel(label + "\n(Å)")
        ax.grid(alpha=0.15)
    axes[0].legend(fontsize=10)
    axes[-1].set_xlabel(
        "Production time (ps); dotted line = dominant grid_3_03 frame at 7.60 ps"
    )
    fig.suptitle(
        "Retained 5–10 ps grid histories: local motion around the dominant frame"
    )
    save(fig, folder, "early_grid_key_windows_5to10ps")
