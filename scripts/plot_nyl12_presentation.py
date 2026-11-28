"""Plot the retained relay PES at presentation scale without changing science."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


def plot_relay(source: Path, destination: Path) -> None:
    """Plot physical energies relative to the regenerated E–S endpoint.

    Args:
        source: Authoritative connected-path CSV.
        destination: Presentation PNG output path.
    """
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    image = [int(row["image"]) for row in rows]
    energy = [float(row["relative_physical_kcal_mol"]) for row in rows]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 29})
    fig, axis = plt.subplots(figsize=(13, 5.0), layout="constrained")
    axis.plot(image, energy, color="#008C95", linewidth=4, marker="o", markersize=5)
    axis.axhline(0, color="#8997A0", linewidth=1)
    axis.set_xlabel("Ordered path configuration", fontsize=29)
    axis.set_ylabel("ΔE (kcal/mol)", fontsize=30)
    axis.set_ylim(-7, 21)
    axis.set_yticks([-5, 0, 5, 10, 15, 20])
    axis.tick_params(labelsize=27)
    axis.spines[["top", "right"]].set_visible(False)
    peak = max(range(len(energy)), key=energy.__getitem__)
    axis.annotate(
        f"{energy[peak]:.2f}",
        (image[peak], energy[peak]),
        xytext=(0, 12),
        textcoords="offset points",
        ha="center",
        fontsize=30,
    )
    axis.annotate(
        "E–S",
        xy=(image[0], energy[0]),
        xytext=(image[-1] * 0.08, -3.8),
        fontsize=29,
        arrowprops={"arrowstyle": "-", "color": "#526577", "linewidth": 1.5},
    )
    axis.annotate(
        "Acylenzyme",
        xy=(image[-1], energy[-1]),
        xytext=(image[-1] * 0.60, -3.8),
        fontsize=29,
        arrowprops={"arrowstyle": "-", "color": "#526577", "linewidth": 1.5},
    )
    fig.savefig(destination, dpi=220, facecolor="white")
    fig.savefig(destination.with_suffix(".svg"), facecolor="white")
    fig.savefig(destination.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    logging.info("Plotted %d source configurations to %s", len(rows), destination)


def plot_wt_occupancy(source: Path, destination: Path) -> None:
    """Plot replica means and standard deviations for matched WT states.

    Args:
        source: Authoritative occupancy summary CSV.
        destination: Presentation PNG output path.
    """
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = {row["state"]: row for row in rows if row["observable"] == "tier_a_or_b"}
    states = ["nh2_asp272_asp274", "nh2_asp272_ash274", "nh3_asp272_asp274"]
    labels = [
        "N-term NH₂\nAsp272 COO⁻\nAsp274 COO⁻",
        "N-term NH₂\nAsp272 COO⁻\nAsp274 COOH",
        "N-term NH₃⁺\nAsp272 COO⁻\nAsp274 COO⁻",
    ]
    means = [float(selected[state]["mean"]) * 100 for state in states]
    deviations = [float(selected[state]["sample_sd"]) * 100 for state in states]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 22})
    figure, axis = plt.subplots(figsize=(12, 5.8), layout="constrained")
    axis.bar(range(3), means, color=["#008C95", "#718A9C", "#17324D"], width=0.64)
    axis.errorbar(
        range(3),
        means,
        yerr=deviations,
        fmt="none",
        ecolor="#223548",
        capsize=8,
        linewidth=2,
    )
    axis.set_xticks(range(3), labels, fontsize=19)
    axis.set_ylabel("Frames meeting geometry criteria (%)", fontsize=23)
    axis.set_ylim(0, 62)
    axis.tick_params(axis="y", labelsize=20)
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#DCE4E9", linewidth=1)
    axis.set_axisbelow(True)
    for index, (mean, deviation) in enumerate(zip(means, deviations, strict=True)):
        axis.text(
            index,
            mean + 2.5,
            f"{mean:.2f} ± {deviation:.2f}%",
            ha="center",
            va="bottom",
            fontsize=18,
        )
    figure.savefig(destination, dpi=220, facecolor="white")
    figure.savefig(destination.with_suffix(".svg"), facecolor="white")
    figure.savefig(destination.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def plot_direct_neb(source: Path, destination: Path) -> None:
    """Plot the canonical Mechanism A NEB energies at presentation scale.

    Args:
        source: Retained NEB image-energy CSV.
        destination: Presentation PNG output path.
    """
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    images = [int(row["image"]) for row in rows]
    energies = [float(row["relative_energy_kcal_mol"]) for row in rows]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 20})
    figure, axis = plt.subplots(figsize=(12, 5.3), layout="constrained")
    axis.plot(images, energies, color="#008C95", linewidth=4, marker="o", markersize=6)
    axis.set_xlabel("NEB image", fontsize=21)
    axis.set_ylabel("Relative potential energy (kcal/mol)", fontsize=21)
    axis.tick_params(labelsize=18)
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#DCE4E9", linewidth=1)
    axis.set_axisbelow(True)
    figure.savefig(destination, dpi=220, facecolor="white")
    figure.savefig(destination.with_suffix(".svg"), facecolor="white")
    figure.savefig(destination.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def plot_direct_pmf(source: Path, destination: Path) -> None:
    """Plot the retained direct-seeded q_acyl marginal and bootstrap interval.

    Args:
        source: Canonical PMF export CSV with bootstrap bounds.
        destination: Presentation PNG output path.
    """
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    coordinate = [float(row["q_acyl_A"]) for row in rows]
    free_energy = [float(row["F_relative_ES_bin_kcal_mol"]) for row in rows]
    lower = [float(row["bootstrap95_low"]) for row in rows]
    upper = [float(row["bootstrap95_high"]) for row in rows]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 20})
    figure, axis = plt.subplots(figsize=(12, 5.3), layout="constrained")
    axis.plot(coordinate, free_energy, color="#008C95", linewidth=4)
    axis.fill_between(
        coordinate,
        lower,
        upper,
        color="#008C95",
        alpha=0.2,
        label="Bootstrap 95% interval",
    )
    axis.set_xlabel(r"Acyl coordinate, $q_{\mathrm{acyl}}$ (Å)", fontsize=21)
    axis.set_ylabel(r"$\Delta G$ (kcal mol$^{-1}$)", fontsize=20)
    axis.tick_params(labelsize=18)
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(color="#DCE4E9", linewidth=1)
    axis.set_axisbelow(True)
    axis.legend(frameon=False, fontsize=16)
    figure.savefig(destination, dpi=220, facecolor="white")
    figure.savefig(destination.with_suffix(".svg"), facecolor="white")
    figure.savefig(destination.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def crop_orientation_panel(source: Path, destination: Path) -> None:
    """Crop old headers and outer whitespace from the paired-orientation panel.

    Args:
        source: Source-bound PyMOL side-by-side orientation figure.
        destination: Curated presentation PNG output path.
    """
    with Image.open(source) as figure:
        width, height = figure.size
        crop = figure.crop(
            (
                int(width * 0.06),
                int(height * 0.12),
                int(width * 0.94),
                int(height * 0.92),
            )
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        crop.save(destination, format="PNG", optimize=True)


def crop_bonded_profiles(source: Path, destination: Path) -> None:
    """Enlarge two representative PA6 bonded-profile panels for the main deck.

    Args:
        source: Eight-panel force-field validation figure.
        destination: Curated two-panel presentation PNG.
    """
    with Image.open(source) as figure:
        width, height = figure.size
        panel_width = width / 4
        top, bottom = int(height * 0.065), int(height * 0.52)
        right_trim = int(panel_width * 0.04)
        crops = [
            figure.crop(
                (
                    int(index * panel_width),
                    top,
                    int((index + 1) * panel_width) - right_trim,
                    bottom,
                )
            )
            for index in (0, 2)
        ]
    panel_width_px, panel_height_px = crops[0].size
    gutter = int(panel_width * 0.16)
    canvas = Image.new(
        "RGB", (panel_width_px * 2 + gutter, panel_height_px), "white"
    )
    canvas.paste(crops[0], (0, 0))
    canvas.paste(crops[1], (panel_width_px + gutter, 0))
    destination.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(destination, format="PNG", optimize=True)


def main() -> None:
    """Render the source-bound relay plot in the existing presentation folder."""
    logging.basicConfig(level=logging.INFO)
    root = Path(__file__).resolve().parents[1]
    plot_relay(
        root
        / "qmmm/pa6_flipped_model30/12_nterminal_relay/plots/connected_relay_pes.csv",
        root / "docs/presentations/nyl12_progress_20261123/assets/relay_pes_large.png",
    )
    plot_wt_occupancy(
        root
        / "qmmm/pa6_flipped_model30/10_wt_protonation_md/analysis_after_md/tables/occupancy_by_state.csv",
        root
        / "docs/presentations/nyl12_progress_20261123/assets/wt_protonation_occupancy.png",
    )
    direct_root = (
        root
        / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113"
    )
    plot_direct_neb(
        root / "docs/neb_us_preliminary_20261121_neb.csv",
        root
        / "docs/presentations/nyl12_progress_20261123/assets/mechanism_A_direct_neb_large.png",
    )
    plot_direct_pmf(
        direct_root
        / "umbrella_sampling/analysis/pmf_refresh_20261122/primary_q_acyl_pmf.csv",
        root
        / "docs/presentations/nyl12_progress_20261123/assets/mechanism_A_direct_seeded_pmf.png",
    )
    crop_orientation_panel(
        root / "analysis/wt_pose_orientation_examples_20261123/side_by_side.png",
        root
        / "docs/presentations/nyl12_progress_20261123/assets/substrate_orientation_comparison_cropped.png",
    )
    crop_bonded_profiles(
        root
        / "docs/presentations/nyl12_progress_20261123/assets/pa6_l4_bonded_profiles.png",
        root
        / "docs/presentations/nyl12_progress_20261123/assets/pa6_l4_bonded_representatives.png",
    )


if __name__ == "__main__":
    main()
