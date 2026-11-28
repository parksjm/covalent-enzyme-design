"""Recover plotting of the retained, converged pooled MBAR fit without refitting."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LOGGER = logging.getLogger(__name__)
WORK = Path(__file__).parent.resolve()


def load_json(path: Path) -> dict[str, Any]:
    """Read a retained JSON evidence record."""
    return json.loads(path.read_text())


def save_json(path: Path, data: dict[str, Any]) -> None:
    """Write finite, explicitly derived JSON evidence."""
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def plot_latest(folder: Path, selection: dict[str, Any]) -> None:
    """Plot unchanged histograms with correctly scoped statistical intervals."""
    primary = load_json(Path(selection["primary_report"]))
    report = load_json(folder / "point.json")
    main = primary["primary"]
    centers = np.asarray(main["centers"], float)
    original = np.asarray(main["pmf_kcal_mol"], float)
    pooled = np.asarray(report["projections"]["hard"]["pmf"], float)
    reference = np.flatnonzero(
        np.isclose(centers, main["reference_s"], rtol=0, atol=1e-12)
    )
    if centers.shape != pooled.shape or len(reference) != 1:
        raise ValueError("Unexpected bins or E–S reference; no remapping allowed")
    if not np.isclose(pooled[reference[0]], 0, atol=1e-10) or not np.isclose(
        original[reference[0]], 0, atol=1e-10
    ):
        raise ValueError("Main and pooled PMFs have different energy references")
    ci = primary["bootstrap"]["5ps"]
    supported = np.asarray(ci["support_fraction"]) == 1
    lower = np.asarray(ci["lower_95_kcal_mol"], float)
    upper = np.asarray(ci["upper_95_kcal_mol"], float)
    lower[~supported] = np.nan
    upper[~supported] = np.nan
    pooled_ess = np.asarray(report["projections"]["hard"]["ESS"], float)
    main_ess = np.asarray(main["weight_ess_per_bin"], float)
    plt.rcParams.update(
        {
            "font.size": 16,
            "axes.labelsize": 18,
            "legend.fontsize": 13,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    fig, (ax, support_ax) = plt.subplots(
        2,
        1,
        figsize=(12, 9),
        sharex=True,
        constrained_layout=True,
        height_ratios=(3.6, 1),
    )
    horizon = selection["main_horizon_ps"]
    ax.fill_between(centers, lower, upper, color="0.5", alpha=0.22)
    ax.plot(
        centers,
        original,
        color="0.35",
        lw=2,
        label=f"Original 57 windows: {horizon} ps + conditional 95% interval",
    )
    ax.plot(
        centers,
        pooled,
        color="tab:orange",
        lw=2.8,
        marker="o",
        ms=4,
        label="All 100 histories: pooled estimate (interval unavailable)",
    )
    ax.axhline(0, color="0.6", lw=1, ls="--")
    ax.scatter([main["reference_s"]], [0], color="0.25", marker="D", s=45, zorder=4)
    ax.set(ylabel=r"PMF relative to E–S (kcal mol$^{-1}$)")
    ax.legend(loc="upper left", frameon=False)
    peak = int(np.nanargmax(pooled))
    ax.annotate(
        f"Projected maximum\n{pooled[peak]:.1f} kcal mol$^{{-1}}$",
        xy=(centers[peak], pooled[peak]),
        xytext=(0.43, 0.44),
        textcoords="axes fraction",
        fontsize=15,
        arrowprops={"arrowstyle": "->", "color": "0.3"},
        bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.95},
    )
    ax.text(
        0.03,
        0.04,
        "Sampling convergence remains unresolved",
        transform=ax.transAxes,
        fontsize=15,
        bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.95},
    )
    support_ax.semilogy(centers, main_ess, color="0.35", lw=1.7)
    support_ax.semilogy(centers, pooled_ess, color="tab:orange", lw=2)
    support_ax.set(
        xlabel=r"Reaction-path progress, $s$",
        ylabel="Weight ESS",
        xlim=(main["edges"][0], main["edges"][-1]),
    )
    support_ax.grid(axis="y", alpha=0.2)
    ax.grid(axis="y", alpha=0.18)
    fig.suptitle("Relay acylation PMF at 300 K", fontsize=23)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(folder / f"relay_pmf_latest.{extension}", dpi=220)
    plt.close(fig)
    fields = [
        "s",
        "pooled_pmf_kcal_mol",
        "main_only_pmf_kcal_mol",
        "main_only_lower95",
        "main_only_upper95",
        "pooled_weight_ess",
    ]
    with (folder / "relay_pmf_latest.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(zip(centers, pooled, original, lower, upper, pooled_ess))
    cohorts: dict[str, Any] = {}
    for lineage in selection["supplemental_lineages"]:
        if lineage["eligible"]:
            cohorts.setdefault(lineage["group"], Counter())[lineage["horizon_ps"]] += 1
    save_json(
        folder / "refresh_summary.json",
        {
            "main_horizon_ps": horizon,
            "main_histories": 57,
            "supplemental_histories_by_horizon_ps": {
                name: dict(counts) for name, counts in cohorts.items()
            },
            "histories": report["states"],
            "retained_frames": report["frames"],
            "pooled_peak_kcal_mol": float(pooled[peak]),
            "peak_s": float(centers[peak]),
            "pooled_CI": "unavailable: one retained 5 ps block/new grid history",
            "gray_interval_scope": "original57 only; conditional pointwise95%; 200 full-MBAR 5ps block resamples",
            "weight_ESS_scope": "bin-local Kish weight concentration, not temporal independence",
            "burn_in": "discard first5ps once per continuous production history",
            "sampling_converged": False,
            "global_barrier_accepted": False,
            "source_publication": str(folder),
            "math_changed": False,
            "curve_smoothing": False,
        },
    )


def main() -> None:
    """Validate and render retained numerical output; never rerun MBAR."""
    logging.basicConfig(level=logging.INFO)
    folder = (WORK / "current").resolve()
    if load_json(folder / "acceptance.json")["process_complete"] is not True:
        raise ValueError("Accepted numerical MBAR fit is required")
    plot_latest(folder, load_json(WORK / "selection.json"))
    payloads = sorted(
        p for p in folder.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt"
    )
    (folder / "SHA256SUMS.txt").write_text(
        "".join(
            f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n"
            for p in payloads
        )
    )
    save_json(
        WORK / "plot_recovery_acceptance.json",
        {
            "process_complete": True,
            "MBAR_refitted": False,
            "original_plot_failure": "E-S reference assumed to be first histogram bin",
            "correct_reference_s": 0.0125,
            "reference_bin_index": 10,
            "all_histogram_bins_displayed": True,
            "math_changed": False,
            "figure": str(folder / "relay_pmf_latest.png"),
        },
    )
    LOGGER.info("Recovered plot from retained MBAR output: %s", folder)


if __name__ == "__main__":
    main()
