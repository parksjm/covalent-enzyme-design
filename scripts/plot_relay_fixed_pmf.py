#!/usr/bin/env python3
"""Render fixed-burn-in PMFs, correlated intervals and analysis sensitivities."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
import logging
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nyl12_covalent.relay_umbrella_analysis_io import digest


def plot_report(report_path: Path, output: Path) -> dict:
    """Plot saved MBAR estimates without modifying scientific selections.

    Args:
        report_path: Fixed-burn-in analysis and bootstrap report.
        output: Campaign plots directory.

    Returns:
        Plot identities and the source report hash.
    """
    data = json.loads(report_path.read_text())
    main, bootstrap = data["primary"], data["bootstrap"]["5ps"]
    x, y = np.array(main["centers"]), np.array(main["pmf_kcal_mol"], float)
    lower = np.array(bootstrap["lower_95_kcal_mol"], float)
    upper = np.array(bootstrap["upper_95_kcal_mol"], float)
    horizon = data["horizon_ps"]
    date = datetime.now(ZoneInfo("America/New_York")).strftime("%Y%m%d")
    prefix = f"relay_{horizon:02.0f}ps_fixed_burnin_{date}"
    output.mkdir(parents=True, exist_ok=True)
    paths = []

    def save(figure: plt.Figure, suffix: str) -> None:
        for extension in ("png", "pdf"):
            path = output / f"{prefix}_{suffix}.{extension}"
            if path.exists():
                raise FileExistsError("Preserve prior figure: " + str(path))
            figure.savefig(path, dpi=200)
            paths.append(str(path))
        plt.close(figure)

    style = {
        "font.size": 18,
        "axes.labelsize": 23,
        "axes.titlesize": 23,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 14,
    }
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(12, 7.5), constrained_layout=True)
        ax.plot(
            x, y, "o-", color="#2366a8", ms=5, lw=2.4, label="All post-burn-in frames"
        )
        ax.fill_between(
            x,
            lower,
            upper,
            color="#2366a8",
            alpha=0.23,
            label="Pointwise 95% CI · 5-ps block bootstrap",
        )
        unsupported = np.isfinite(y) & ~np.isfinite(lower)
        if unsupported.any():
            ax.plot(
                x[unsupported],
                y[unsupported],
                "o",
                color="0.4",
                ms=5,
                label="CI unresolved: incomplete bootstrap support",
            )
        ax.plot(main["reference_s"], 0, "D", color="#d17b00", ms=8)
        ax.axhline(0, color="0.6", lw=1, ls="--")
        ax.set(
            xlabel="Normalized path coordinate, $s$",
            ylabel="PMF relative to E–S reference\n(kcal mol⁻¹)",
            title=f"N-terminal relay · {horizon:g} ps/window · 300 K",
        )
        ax.text(
            0.02,
            0.97,
            "5-ps production burn-in\nSampling not converged",
            transform=ax.transAxes,
            va="top",
            fontsize=15,
            bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.95),
        )
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="lower left")
        save(fig, "pmf")

        fig, axes = plt.subplots(
            1, 2, figsize=(16, 7), constrained_layout=True, sharey=True
        )
        axes[0].plot(x, y, color="#2366a8", lw=2.8, label="5-ps burn-in (primary)")
        for key, color in zip(
            ("burn_0ps", "burn_2.5ps", "burn_10ps"),
            ("#999999", "#009e73", "#d55e00"),
            strict=True,
        ):
            item = data["sensitivity"][key]
            axes[0].plot(
                x,
                np.array(item["pmf_kcal_mol"], float),
                color=color,
                lw=1.7,
                label=f"{item['burn_in_ps']:g}-ps burn-in",
            )
        axes[0].set_title("Burn-in sensitivity")
        midpoint = (5 + horizon) / 2
        for key, color, label in (
            ("early_half", "#009e73", f"5–{midpoint:g} ps"),
            ("late_half", "#d55e00", f"{midpoint:g}–{horizon:g} ps"),
        ):
            axes[1].plot(
                x,
                np.array(data["sensitivity"][key]["pmf_kcal_mol"], float),
                color=color,
                lw=2,
                label=label,
            )
        axes[1].plot(
            x, y, color="#2366a8", lw=2, ls="--", label="Combined post-burn-in"
        )
        axes[1].set_title("Early versus late sampling")
        for ax in axes:
            ax.set_xlabel("Normalized path coordinate, $s$")
            ax.grid(axis="y", alpha=0.2)
            ax.legend(loc="lower left")
        axes[0].set_ylabel("PMF relative to fixed E–S bin\n(kcal mol⁻¹)")
        save(fig, "sensitivity")

        fig, ax = plt.subplots(figsize=(12, 7.5), constrained_layout=True)
        ax.plot(x, y, color="#2366a8", lw=2)
        for key, color, alpha in (("5ps", "#2366a8", 0.22), ("2.5ps", "#d55e00", 0.15)):
            item = data["bootstrap"][key]
            ax.fill_between(
                x,
                np.array(item["lower_95_kcal_mol"], float),
                np.array(item["upper_95_kcal_mol"], float),
                color=color,
                alpha=alpha,
                label=f"{key} blocks · pointwise 95% CI",
            )
        ax.set(
            xlabel="Normalized path coordinate, $s$",
            ylabel="PMF relative to E–S reference\n(kcal mol⁻¹)",
            title=f"Uncertainty depends on block length · {horizon:g} ps/window",
        )
        ax.grid(axis="y", alpha=0.2)
        ax.legend()
        save(fig, "block_sensitivity")
    metadata = dict(
        source=str(report_path.resolve()),
        source_sha256=digest(report_path),
        plot_source_sha256=digest(Path(__file__)),
        plots=paths,
        confidence_method=bootstrap["method"],
        no_smoothing=True,
    )
    (output / f"{prefix}_source.json").write_text(json.dumps(metadata, indent=2) + "\n")
    logging.info("Rendered fixed-burn-in PMF and sensitivities: %s", paths[0])
    return metadata


def main() -> None:
    """Read saved numerical evidence and render consistent PNG/PDF figures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    plot_report(args.report, args.output_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
