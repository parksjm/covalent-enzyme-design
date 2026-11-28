"""Render exact movie distances and native biased-window distributions."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from nyl12_covalent.reaction_distance_stats import (  # noqa: E402
    DISTANCES,
    LABELS,
    ROLES,
    campaign_metadata,
    load_sources,
    movie_rows,
    ownership_summary,
    sha256,
    summarize,
    write_csv,
)

LOGGER = logging.getLogger(__name__)


def render(
    movie: list[dict],
    rows: list[dict],
    windows: list[str],
    out: Path,
    charts: dict[str, str],
) -> None:
    """Draw independent-window traces and descriptive distribution panels."""
    meta = campaign_metadata(movie, rows, windows, charts)
    groups = meta["windows"]

    def span(key: str) -> str:
        values = [g[key] for g in groups]
        return (
            f"{min(values):g}"
            if min(values) == max(values)
            else f"{min(values):g}–{max(values):g}"
        )

    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(
        5, 1, figsize=(17, 14), sharex=True, constrained_layout=True
    )
    for key, label, ax in zip(DISTANCES, LABELS, axes, strict=True):
        for wi, window in enumerate(windows):
            selected = [r for r in movie if r["window"] == window]
            x = np.array([r["movie_frame0based"] for r in selected])
            y = np.array([r[key] for r in selected])
            ax.scatter(
                x, y, s=5, color="#0072B2", alpha=0.55, linewidths=0, rasterized=True
            )
            # Within-window running median; independent windows never join.
            median = [
                np.median(y[max(0, i - 2) : min(len(y), i + 3)]) for i in range(len(y))
            ]
            ax.plot(x, median, color="black", lw=0.65)
            ax.axvline(x[0] - 0.5, color=".85", lw=0.4)
        for boundary in meta["phase_boundary_movie_frames"]:
            ax.axvline(boundary, color="#D55E00", lw=1.5, linestyle="--")
        ax.set_ylabel("Distance (Å)")
        ax.set_title(label, loc="left", fontweight="bold", fontsize=14)
        ax.grid(axis="y", alpha=0.18)
    for chart in dict.fromkeys(charts[w] for w in windows):
        selected = [g for g in groups if g["chart"] == chart]
        center = (selected[0]["movie_start"] + selected[-1]["movie_stop"]) / (
            2 * len(movie)
        )
        axes[0].text(
            center,
            1.21,
            f"{chart.title()} · {len(selected)} independent windows",
            transform=axes[0].transAxes,
            ha="center",
        )
    axes[-1].set_xlabel(
        f"Path-ordered movie frame (0-based; {span('movie_count')} frames/window; not MD time)"
    )
    fig.suptitle(
        f"Relay B · exact {span('horizon_ps')} ps/window VMD movie distances\n"
        "Blue: stored movie frames · Black: 5-frame within-window median",
        fontsize=18,
    )
    for extension in ("png", "pdf"):
        fig.savefig(out / f"reaction_distances.{extension}", dpi=220)
    plt.close(fig)
    fig, axes = plt.subplots(
        5, 1, figsize=(17, 14), sharex=True, constrained_layout=True
    )
    for key, label, ax in zip(DISTANCES, LABELS, axes, strict=True):
        values = [
            np.array(
                [
                    r[key]
                    for r in rows
                    if r["window"] == w and r["production_time_ps"] > 5
                ]
            )
            for w in windows
        ]
        parts = ax.violinplot(
            values,
            positions=np.arange(1, len(windows) + 1),
            widths=0.85,
            showextrema=False,
            showmedians=True,
            points=60,
        )
        for body in parts["bodies"]:
            body.set_facecolor("#0072B2")
            body.set_alpha(0.6)
        parts["cmedians"].set_color("black")
        for boundary in meta["phase_boundary_window_positions"]:
            ax.axvline(boundary, color="#D55E00", linestyle="--", lw=1.5)
        ax.set(title=label, ylabel="Distance (Å)")
        ax.grid(axis="y", alpha=0.18)
    stride = max(1, int(np.ceil(len(windows) / 12)))
    ticks = sorted(set(range(1, len(windows) + 1, stride)) | {len(windows)})
    axes[-1].set_xticks(
        ticks,
        [
            f"{charts[windows[i - 1]][0].upper()}{windows[i - 1].rsplit('_', 1)[-1]}"
            for i in ticks
        ],
    )
    axes[-1].set_xlabel("Independent biased window in path order")
    fig.suptitle(
        f"Native distances · {span('retained_count')} frames/window after one 5 ps discard\n"
        "Descriptive biased distributions · black bars: medians",
        fontsize=18,
    )
    for extension in ("png", "pdf"):
        fig.savefig(out / f"per_window_distributions.{extension}", dpi=220)
    plt.close(fig)


def main() -> None:
    """Build reproducible distance plots and tabular descriptive statistics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    movie_dir = args.base / "vmd_progress/current"
    manifest = json.loads((movie_dir / "manifest.json").read_text())
    if any(manifest["roles_1based"][k] != v for k, v in ROLES.items()):
        raise ValueError("Unexpected original atom roles")
    args.output.mkdir(parents=True, exist_ok=False)
    rows, sources = load_sources(args.base, manifest)
    LOGGER.info("Read %d native frames", len(rows))
    movie = movie_rows(rows, manifest["frame_plan"])
    stats = summarize(rows, manifest["windows"])
    meta = campaign_metadata(
        movie, rows, manifest["windows"], manifest["window_charts"]
    )
    write_csv(args.output / "movie_frame_distances.csv", movie)
    write_csv(
        args.output / "native_retained_distances.csv",
        [
            {
                "window": r["window"],
                "production_time_ps": r["production_time_ps"],
                **{key: r[key] for key in DISTANCES},
            }
            for r in rows
            if r["production_time_ps"] > 5
        ],
    )
    write_csv(args.output / "per_window_stats.csv", stats)
    ownership = ownership_summary(rows, manifest["windows"])
    (args.output / "per_window_stats.json").write_text(
        json.dumps({"statistics": stats, "hydrogen_ownership": ownership}, indent=2)
        + "\n"
    )
    sources[str((movie_dir / "manifest.json").resolve())] = sha256(
        movie_dir / "manifest.json"
    )
    for name in (
        "src/nyl12_covalent/reaction_distance_stats.py",
        "scripts/plot_relay_reaction_distance_statistics.py",
    ):
        sources[str(Path(name).resolve())] = sha256(Path(name))
    (args.output / "source_provenance.json").write_text(
        json.dumps(
            {
                "sha256": sources,
                "roles_original_1based": ROLES,
                "movie_manifest_trajectory_sha256": manifest["trajectory_sha256"],
                "discard_ps_once": 5,
                **meta,
            },
            indent=2,
        )
        + "\n"
    )
    render(movie, rows, manifest["windows"], args.output, manifest["window_charts"])
    LOGGER.info("Wrote plots and descriptive statistics to %s", args.output)


if __name__ == "__main__":
    main()
