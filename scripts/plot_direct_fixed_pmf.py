#!/usr/bin/env python3
"""Plot saved fixed-window direct and relay PMF reports without refitting."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
import logging
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _digest(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _values(values: list[float | None]) -> np.ndarray:
    """Convert nullable report values to floats while preserving null gaps."""
    return np.asarray([np.nan if value is None else value for value in values], float)


def _axes(
    ax: plt.Axes, route: str, horizon: float, reference: float, temperature: float
) -> None:
    """Apply consistent labels and route-specific arclength explanation."""
    ax.plot(reference, 0, "D", color="#d17b00", ms=8)
    ax.axhline(0, color="0.6", lw=1, ls="--")
    ax.set(
        xlabel="Route-specific normalized path coordinate, $s$",
        ylabel="PMF relative to E–S\n(kcal mol⁻¹)",
        title=f"{route}\n{horizon:g} ps/window · {temperature:g} K",
    )
    ax.grid(axis="y", alpha=0.2)


def _save_pair(fig: plt.Figure, output: Path, stem: str) -> list[str]:
    """Save a figure as new PNG and PDF artifacts."""
    paths = [output / f"{stem}.{ext}" for ext in ("png", "pdf")]
    existing = [path for path in paths if path.exists()]
    if existing:
        raise FileExistsError(
            "Preserve prior figures: " + ", ".join(map(str, existing))
        )
    for path in paths:
        fig.savefig(path, dpi=200)
    plt.close(fig)
    return [str(path) for path in paths]


def _draw_main(ax: plt.Axes, report: dict, route: str, color: str) -> None:
    """Draw a report's primary PMF and pointwise confidence band."""
    main = report["primary"]
    boot = report["bootstrap"]["5ps"]
    x = np.asarray(main["centers"], float)
    y = _values(main["pmf_kcal_mol"])
    lower = _values(boot["lower_95_kcal_mol"])
    upper = _values(boot["upper_95_kcal_mol"])
    ax.plot(x, y, "o-", color=color, ms=5, lw=2.4, label="MBAR estimate")
    ax.fill_between(
        x,
        lower,
        upper,
        color=color,
        alpha=0.23,
        label="Pointwise 95% CI · 5-ps blocks",
    )
    unsupported = np.isfinite(y) & ~np.isfinite(lower)
    if unsupported.any():
        ax.plot(
            x[unsupported],
            y[unsupported],
            "o",
            color="0.4",
            ms=5,
            label="CI unresolved",
        )
    _axes(
        ax,
        route,
        float(report["horizon_ps"]),
        float(main["reference_s"]),
        float(report["temperature_K"]),
    )
    ax.legend(loc="upper right")


def _draw_sensitivity(axes: np.ndarray, report: dict) -> None:
    """Draw requested burn-in and early/late sensitivity estimates."""
    main = report["primary"]
    data = report["sensitivity"]
    horizon = float(report["horizon_ps"])
    x = np.asarray(main["centers"], float)
    axes[0].plot(
        x,
        _values(main["pmf_kcal_mol"]),
        color="#2366a8",
        lw=2.8,
        label="5-ps burn-in (primary)",
    )
    for key, color in (
        ("burn_0ps", "#999999"),
        ("burn_2.5ps", "#009e73"),
        ("burn_10ps", "#d55e00"),
    ):
        item = data[key]
        axes[0].plot(
            x,
            _values(item["pmf_kcal_mol"]),
            color=color,
            lw=1.7,
            label=f"{item['burn_in_ps']:g}-ps burn-in",
        )
    axes[0].set_title("Burn-in sensitivity")
    midpoint = (5.0 + horizon) / 2
    for key, color, label in (
        ("early_half", "#009e73", f"5–{midpoint:g} ps"),
        ("late_half", "#d55e00", f"{midpoint:g}–{horizon:g} ps"),
    ):
        axes[1].plot(
            x, _values(data[key]["pmf_kcal_mol"]), color=color, lw=2, label=label
        )
    axes[1].plot(
        x,
        _values(main["pmf_kcal_mol"]),
        color="#2366a8",
        lw=2,
        ls="--",
        label="Combined post-burn-in",
    )
    axes[1].set_title("Early versus late sampling")
    for ax in axes:
        ax.set(
            xlabel="Route-specific normalized path coordinate, $s$",
            ylabel="PMF relative to E–S\n(kcal mol⁻¹)",
        )
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="lower left")


def plot_reports(direct_path: Path, relay_path: Path, output: Path) -> dict:
    """Render report-only direct PMFs, sensitivities, and route comparison.

    Args:
        direct_path: Saved direct-route fixed-PMF report.
        relay_path: Saved relay-route fixed-PMF report.
        output: Destination directory for new figures and metadata.

    Returns:
        Artifact paths, input hashes, and the matched-horizon designation.
    """
    direct = json.loads(direct_path.read_text())
    relay = json.loads(relay_path.read_text())
    if direct["protocol"] != relay["protocol"]:
        raise ValueError("Comparison requires the same analysis protocol")
    if direct["temperature_K"] != relay["temperature_K"]:
        raise ValueError("Comparison requires the same physical temperature")
    direct_h, relay_h = float(direct["horizon_ps"]), float(relay["horizon_ps"])
    matched = np.isclose(direct_h, relay_h)
    date = datetime.now(ZoneInfo("America/New_York")).strftime("%Y%m%d")
    direct_stem = f"direct_{direct_h:02.0f}ps_fixed_{date}"
    comparison_stem = f"comparison_{direct_h:g}ps_vs_relay_{relay_h:g}ps_{date}"
    output.mkdir(parents=True, exist_ok=True)
    style = {
        "font.size": 18,
        "axes.labelsize": 21,
        "axes.titlesize": 22,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "legend.fontsize": 13,
    }
    artifacts: list[str] = []
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(12, 7.5), constrained_layout=True)
        _draw_main(ax, direct, "Direct-NEB-seeded sampling", "#2366a8")
        ax.text(
            0.02,
            0.97,
            "Proton channels not restricted\n5-ps production burn-in\n"
            "Sampling not converged",
            transform=ax.transAxes,
            va="top",
            fontsize=15,
            bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.95),
        )
        artifacts += _save_pair(fig, output, direct_stem + "_pmf")

        fig, axes = plt.subplots(
            1, 2, figsize=(16, 7), constrained_layout=True, sharey=True
        )
        _draw_sensitivity(axes, direct)
        fig.suptitle(
            "Direct-NEB-seeded sampling · sensitivity · sampling not converged"
        )
        artifacts += _save_pair(fig, output, direct_stem + "_sensitivity")

        fig, axes = plt.subplots(
            1, 2, figsize=(17, 7.5), constrained_layout=True, sharey=True
        )
        _draw_main(axes[0], direct, "Direct-NEB-seeded sampling", "#2366a8")
        _draw_main(axes[1], relay, "N-terminal-relay-seeded sampling", "#d55e00")
        if matched and np.isclose(direct_h, 25):
            comparison_title = "Same analysis protocol; sampling not converged"
        elif matched:
            comparison_title = (
                f"Matched {direct_h:g}-ps horizons; sampling not converged"
            )
        else:
            comparison_title = (
                f"Unmatched horizons ({direct_h:g} vs {relay_h:g} ps); "
                "sampling not converged"
            )
        fig.suptitle(comparison_title)
        artifacts += _save_pair(fig, output, comparison_stem)

    metadata = {
        "direct_report": str(direct_path.resolve()),
        "direct_report_sha256": _digest(direct_path),
        "relay_report": str(relay_path.resolve()),
        "relay_report_sha256": _digest(relay_path),
        "plot_script": str(Path(__file__).resolve()),
        "plot_script_sha256": _digest(Path(__file__)),
        "direct_horizon_ps": direct_h,
        "relay_horizon_ps": relay_h,
        "horizons_matched": bool(matched),
        "horizon_comparison": "matched25both"
        if matched and direct_h == 25
        else ("matched" if matched else f"unmatched{direct_h:g}vs{relay_h:g}"),
        "arclength_note": (
            "Route-specific normalized arclength uses the same 3-CV metric and "
            "projection, bin width, and fixed reference."
        ),
        "sampling_converged": False,
        "no_smoothing": True,
        "plots": artifacts,
    }
    metadata_path = output / f"{comparison_stem}_source.json"
    if metadata_path.exists():
        raise FileExistsError("Preserve prior metadata: " + str(metadata_path))
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    logging.info("Rendered direct fixed-PMF artifacts: %s", artifacts[0])
    return metadata


def main() -> None:
    """Parse CLI arguments and render figures from saved reports."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct-report", type=Path, required=True)
    parser.add_argument("--relay-report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    plot_reports(args.direct_report, args.relay_report, args.output_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
