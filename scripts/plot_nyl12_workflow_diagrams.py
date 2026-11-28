#!/usr/bin/env python3
"""Render editable mechanism and modeling workflow schematics for the NYL12 deck."""

from __future__ import annotations

import logging
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/presentations/nyl12_progress_20261123/assets"
NAVY = "#12233F"
BLUE = "#DCEAF7"
PALE = "#F3F6FA"
ACCENT = "#2671A8"


def add_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    wh: tuple[float, float],
    label: str,
    *,
    fontsize: float = 25,
    dashed: bool = False,
    fill: str = "white",
) -> None:
    """Add a rounded text box with a clearly visible outline.

    Args:
        ax: Axes receiving the box.
        xy: Lower-left coordinates in axes data space.
        wh: Box width and height.
        label: Text centered inside the box.
        fontsize: Text size in points.
        dashed: Whether to use a dashed outline for a proposed intermediate.
        fill: Face color for the box.
    """
    box = FancyBboxPatch(
        xy,
        *wh,
        boxstyle="round,pad=0.12,rounding_size=0.12",
        linewidth=2.0,
        edgecolor=NAVY,
        facecolor=fill,
        linestyle="--" if dashed else "-",
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + wh[0] / 2,
        xy[1] + wh[1] / 2,
        label,
        ha="center",
        va="center",
        color=NAVY,
        fontsize=fontsize,
        linespacing=1.25,
    )


def add_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    label: str,
    *,
    fontsize: float = 23,
    color: str = NAVY,
) -> None:
    """Draw a labeled straight reaction or workflow arrow.

    Args:
        ax: Axes receiving the arrow.
        start: Arrow tail coordinates.
        end: Arrow head coordinates.
        label: Short process label centered above the arrow.
        fontsize: Label size in points.
        color: Arrow and label color.
    """
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=25,
        linewidth=2.2,
        color=color,
    )
    ax.add_patch(arrow)
    ax.text(
        (start[0] + end[0]) / 2,
        (start[1] + end[1]) / 2 + 0.22,
        label,
        ha="center",
        va="bottom",
        color=color,
        fontsize=fontsize,
    )


def mechanism_figure() -> plt.Figure:
    """Build two independent, left-to-right acylation channel rows."""
    fig, ax = plt.subplots(figsize=(18, 9))
    ax.set(xlim=(0, 18), ylim=(0, 9))
    ax.axis("off")
    fig.patch.set_facecolor("white")
    # Duplicate endpoint labels denote the same chemical states, not one shared pose.
    common_es = "COMMON E–S\nThr233 OG1–H\nNαH₂"
    product = "COMMON\nPRODUCT\nAcyl–Thr233\n+ amine product"

    ax.text(
        0.45, 8.25, "A  Direct proton transfer", fontsize=27, weight="bold", color=NAVY
    )
    add_box(ax, (0.45, 5.2), (3.0, 1.65), common_es, fontsize=24, fill=BLUE)
    add_box(ax, (14.55, 5.2), (3.0, 1.65), product, fontsize=24, fill=BLUE)
    add_arrow(ax, (3.55, 6.0), (14.45, 6.0), "", fontsize=1)
    ax.text(
        9.0,
        6.55,
        "H: OG1 → amide N\nOG1 attack at C; C–N cleavage",
        ha="center",
        va="bottom",
        fontsize=24,
        color=ACCENT,
    )

    ax.text(
        0.45, 4.45, "B  N-terminal proton relay", fontsize=27, weight="bold", color=NAVY
    )
    add_box(ax, (0.45, 1.25), (3.0, 1.65), common_es, fontsize=24, fill=BLUE)
    add_box(
        ax,
        (6.9, 1.25),
        (3.8, 1.65),
        "Addition /\nproton-transfer region\nNα temporarily NH₃⁺",
        fontsize=22,
        dashed=True,
        fill=PALE,
    )
    add_box(ax, (14.55, 1.25), (3.0, 1.65), product, fontsize=24, fill=BLUE)
    add_arrow(ax, (3.55, 2.05), (6.8, 2.05), "", fontsize=1)
    ax.text(
        5.18,
        2.6,
        "H: OG1 → Nα\nC–OG1 formation",
        ha="center",
        va="bottom",
        fontsize=20,
        color=ACCENT,
    )
    add_arrow(ax, (10.8, 2.05), (14.45, 2.05), "", fontsize=1)
    ax.text(
        12.63,
        2.55,
        "H: Nα → amide N\nC–N cleavage",
        ha="center",
        va="bottom",
        fontsize=20,
        color=ACCENT,
    )
    ax.text(
        0.45,
        0.55,
        "Nα = Thr233 N-terminal amine; amide N = leaving-group nitrogen. No stable TI1 minimum is implied.\n"
        "Proton transfer, attack, and cleavage may be asynchronous; acyl-enzyme hydrolysis is not shown.",
        fontsize=16,
        color=NAVY,
    )
    fig.subplots_adjust(left=0.025, right=0.985, top=0.98, bottom=0.07)
    return fig


def workflow_figure() -> plt.Figure:
    """Build an eight-stage serpentine workflow and benchmark crossbar."""
    fig, ax = plt.subplots(figsize=(18, 9))
    ax.set(xlim=(0, 18), ylim=(0, 9))
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.text(
        0.45,
        8.35,
        "Design objective: increased catalytic efficiency",
        fontsize=24,
        color=NAVY,
    )
    top_row = [
        "PA6 RESP /\nGAFF2 params",
        "Classical WT MD /\nprotonation tests",
        "MM + QM/MM\nE–S minimization",
        "Reaction paths\nCoordinate driving\nNEB",
    ]
    bottom_row = [
        "Boltz-2\nOpenMM",
        "LigandMPNN",
        "TS verification /\ntemplates",
        "Free energies\nUmbrella sampling\nMBAR / WHAM",
    ]
    x_positions = [0.45, 4.95, 9.45, 13.95]
    y_top, y_bottom = 6.15, 2.05
    bw, bh = 3.45, 1.55
    for col, x in enumerate(x_positions):
        add_box(ax, (x, y_top), (bw, bh), top_row[col], fontsize=24, fill=BLUE)
        add_box(ax, (x, y_bottom), (bw, bh), bottom_row[col], fontsize=24, fill="white")
        if col < 3:
            add_arrow(
                ax,
                (x + bw + 0.08, y_top + bh / 2),
                (x_positions[col + 1] - 0.12, y_top + bh / 2),
                "",
                fontsize=1,
            )

    # Top row flows left-to-right; bottom row continues from right-to-left.
    add_arrow(ax, (15.68, y_top - 0.12), (15.68, y_bottom + bh + 0.12), "", fontsize=1)
    for col in range(3, 0, -1):
        start_x = x_positions[col] - 0.12
        end_x = x_positions[col - 1] + bw + 0.12
        add_arrow(
            ax, (start_x, y_bottom + bh / 2), (end_x, y_bottom + bh / 2), "", fontsize=1
        )

    add_box(
        ax,
        (6.35, 4.05),
        (6.2, 1.35),
        "QM method benchmarks\nQM-region selection\nInform QM/MM setup",
        fontsize=22,
        dashed=True,
        fill=PALE,
    )
    add_arrow(ax, (10.2, 5.42), (10.2, y_top - 0.12), "", fontsize=1, color=ACCENT)
    add_arrow(ax, (2.18, y_bottom - 0.12), (2.18, 1.0), "", fontsize=1)
    ax.text(
        0.45,
        0.55,
        "Shortlist for physical + experimental validation",
        fontsize=24,
        color=NAVY,
        ha="left",
        va="center",
    )
    fig.subplots_adjust(left=0.025, right=0.985, top=0.98, bottom=0.04)
    return fig


def save_figure(fig: plt.Figure, stem: str) -> None:
    """Save PNG, SVG, and PDF exports for one figure.

    Args:
        fig: Figure to export.
        stem: Output filename stem in the assets directory.
    """
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg", "pdf"):
        path = OUTPUT / f"{stem}.{suffix}"
        fig.savefig(path, dpi=220, facecolor="white")
        LOGGER.info("Wrote %s", path)
    plt.close(fig)


def main() -> None:
    """Render both requested diagrams in PNG and editable vector formats."""
    save_figure(mechanism_figure(), "mechanism_A_B_schematic")
    save_figure(workflow_figure(), "modeling_workflow")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
