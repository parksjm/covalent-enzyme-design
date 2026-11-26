"""Individual report-only relay PMF and overlap plot panels."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

LOGGER = logging.getLogger(__name__)


def finite_mask(values: np.ndarray) -> np.ndarray:
    """Return the finite-value mask used to retain empty histogram bins."""
    return np.isfinite(values)


def save_pair(figure: plt.Figure, destination: Path) -> list[Path]:
    """Save a figure in PNG and PDF formats and close it."""
    paths = [destination.with_suffix(".png"), destination.with_suffix(".pdf")]
    for path in paths:
        figure.savefig(path, dpi=200 if path.suffix == ".png" else None)
    plt.close(figure)
    return paths


def plot_pmf(
    data: dict[str, Any],
    output: Path,
    stem: str,
    report_sha256: str,
    interval_method: str,
) -> list[Path]:
    """Draw the raw PMF, selected interval, and separate support panels."""
    arrays = data["arrays"]
    centers = arrays["centers"]
    pmf = arrays["pmf_kcal_mol"]
    selected = arrays["selected_frames_per_bin"]
    ess = arrays["weight_ess_per_bin"]
    interval_label: str
    if interval_method == "analytic":
        lower = arrays["lower_95_kcal_mol"]
        upper = arrays["upper_95_kcal_mol"]
        interval_label = "Pointwise approximate 95% asymptotic CI"
        interval_note = "Asymptotic normal CI may understate uncertainty at low ESS"
        support_fraction = 1.0
    elif interval_method == "block-bootstrap":
        bootstrap = data["uncertainty"].get("block_bootstrap")
        if bootstrap is None:
            raise ValueError("Block-bootstrap intervals were requested but are absent")
        lower = arrays["block_bootstrap_lower_95_kcal_mol"]
        upper = arrays["block_bootstrap_upper_95_kcal_mol"]
        support_fraction = arrays["block_bootstrap_support_fraction"]
        supported_bins = np.isfinite(support_fraction) & (support_fraction >= 1)
        interval_label = (
            f"Pointwise 95% 5-ps block-bootstrap CI "
            f"({bootstrap['replicates']} MBAR refits)"
        )
        interval_note = "Block-bootstrap CI; conditional statistical uncertainty only"
    else:
        raise ValueError("interval_method must be analytic or block-bootstrap")

    figure, (axis, counts_axis, ess_axis) = plt.subplots(
        3,
        1,
        figsize=(14, 12),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3.2, 1, 1]},
        facecolor="white",
    )
    valid = finite_mask(pmf)
    interval = valid & finite_mask(lower) & finite_mask(upper)
    if interval_method == "block-bootstrap":
        interval &= supported_bins
    axis.plot(
        centers,
        np.where(valid, pmf, np.nan),
        "o-",
        color="#2366a8",
        lw=2.2,
        ms=5,
        label="Unsmoothed MBAR PMF",
    )
    axis.fill_between(
        centers,
        lower,
        upper,
        where=interval,
        interpolate=False,
        color="#7fb3d5",
        alpha=0.35,
        label=interval_label if interval.any() else None,
    )
    axis.plot(
        data["reference_s"],
        0,
        "D",
        color="#d17b00",
        ms=10,
        label="Common E–S reference",
    )
    axis.axhline(0, color="0.5", lw=1, ls="--")
    axis.set_title(
        f"N-terminal relay PMF at {data['horizon_ps']:g} ps per window · "
        f"{data['temperature_K']:g} K · sampling unconverged"
    )
    axis.set_ylabel("PMF relative to common E–S reference\n(kcal mol⁻¹)")
    axis.text(
        0.015,
        0.97,
        f"{interval_note}\n{_block_count_note(data)} cannot exclude unsampled slow states\n"
        "Red bars: local weight ESS < 5; empty PMF bins are not interpolated",
        transform=axis.transAxes,
        va="top",
        fontsize=14,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.92},
    )
    axis.legend(loc="best", fontsize=17)
    edges = np.asarray(data["report"]["pmf"]["edges"], dtype=float)
    width = np.diff(edges) * 0.82
    selected_finite = finite_mask(selected)
    ess_finite = finite_mask(ess)
    counts_axis.bar(
        centers[selected_finite],
        selected[selected_finite],
        width=width[selected_finite],
        color="#668cae",
    )
    low_ess = ess_finite & (ess < 5)
    colors = np.where(low_ess[ess_finite], "#bb4430", "#d89b48")
    ess_axis.bar(
        centers[ess_finite], ess[ess_finite], width=width[ess_finite], color=colors
    )
    ess_axis.axhline(5, color="#bb4430", ls="--", lw=1.2)
    counts_axis.set_ylabel("Selected\nframes / bin")
    ess_axis.set_ylabel("Local weight\nESS / bin")
    ess_axis.set_xlabel("Normalized path coordinate, $s$")
    ess_axis.set_xlim(edges[0], edges[-1])
    for panel in (axis, counts_axis, ess_axis):
        panel.grid(axis="y", alpha=0.2)
        panel.set_axisbelow(True)
    paths = save_pair(figure, output / f"{stem}_pmf")
    LOGGER.info("Rendered PMF figure from analysis report SHA-256 %s", report_sha256)
    return paths


def _block_count_note(data: dict[str, Any]) -> str:
    """Describe the accepted source-block count from this report."""
    counts = data.get("blocks_per_window_range")
    if counts is None:
        return "Source-block count unavailable"
    if counts[0] == counts[1]:
        return f"{counts[0]} accepted source blocks/window"
    return f"{counts[0]}–{counts[1]} accepted source blocks/window"


def plot_overlap(data: dict[str, Any], output: Path, stem: str) -> list[Path]:
    """Draw the full MBAR overlap matrix using the report's ordered states."""
    matrix = data["overlap"]
    count = len(data["labels"])
    figure, axis = plt.subplots(
        figsize=(14, 12), constrained_layout=True, facecolor="white"
    )
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")
    artist = axis.imshow(
        np.ma.masked_equal(matrix, 0),
        cmap=cmap,
        origin="upper",
        norm=LogNorm(vmin=1e-6, vmax=1, clip=True),
        interpolation="nearest",
    )
    stride = max(1, int(np.ceil(count / 12)))
    tick_indices = np.arange(0, count, stride)
    if tick_indices[-1] != count - 1:
        tick_indices = np.append(tick_indices, count - 1)
    axis.set_xticks(
        tick_indices,
        labels=[data["labels"][i] for i in tick_indices],
        rotation=45,
        ha="right",
    )
    axis.set_yticks(tick_indices, labels=[data["labels"][i] for i in tick_indices])
    axis.set(xlabel="Umbrella window j", ylabel="Umbrella window i")
    axis.set_title(
        f"Full {count}-window MBAR overlap at {data['horizon_ps']:g} ps/window · "
        f"{data['temperature_K']:g} K · sampling unconverged"
    )
    bar = figure.colorbar(artist, ax=axis, shrink=0.82, pad=0.025)
    bar.set_label("MBAR overlap $O_{ij}$ (log scale)")
    paths = save_pair(figure, output / f"{stem}_mbar_overlap")
    return paths


def plot_neighbors(data: dict[str, Any], output: Path, stem: str) -> list[Path]:
    """Plot sorted pair intersections without connecting points as a curve."""
    pairs = data["pairs"]
    values = data["intersections"]
    midpoints = pairs.mean(axis=1) + 1
    order = np.argsort(midpoints, kind="stable")
    pairs, values, midpoints = pairs[order], values[order], midpoints[order]
    adjacent = np.abs(pairs[:, 1] - pairs[:, 0]) == 1
    passes = values >= data["overlap_threshold"]
    figure, axis = plt.subplots(
        figsize=(14, 8), constrained_layout=True, facecolor="white"
    )
    styles = (
        (adjacent & passes, "o", "#2366a8", "Adjacent, above threshold"),
        (adjacent & ~passes, "o", "#bb4430", "Adjacent, below threshold"),
        (~adjacent & passes, "D", "#2366a8", "Nonadjacent, above threshold"),
        (~adjacent & ~passes, "D", "#bb4430", "Nonadjacent, below threshold"),
    )
    for mask, marker, color, label in styles:
        if not mask.any():
            continue
        axis.scatter(
            midpoints[mask],
            values[mask],
            marker=marker,
            color=color,
            s=75,
            label=label,
            zorder=3,
        )
        if marker == "D":
            for x, y, pair in zip(
                midpoints[mask], values[mask], pairs[mask], strict=True
            ):
                axis.annotate(
                    f"{pair[0] + 1}–{pair[1] + 1}",
                    (x, y),
                    xytext=(5, 7),
                    textcoords="offset points",
                    fontsize=12,
                )
    axis.axhline(
        data["overlap_threshold"],
        color="#bb4430",
        ls="--",
        lw=1.5,
        label="Screening threshold",
    )
    axis.set(
        xlabel="Pair midpoint (manifest state number)",
        ylabel="3CV histogram intersection",
        title=f"Window overlap screening at {data['horizon_ps']:g} ps/window · {data['temperature_K']:g} K",
        ylim=(-0.035, 1.04),
    )
    axis.grid(axis="y", alpha=0.2)
    axis.legend(fontsize=15, loc="best")
    paths = save_pair(figure, output / f"{stem}_neighbor_overlap")
    return paths
