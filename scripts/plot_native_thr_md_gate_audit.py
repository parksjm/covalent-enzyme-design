"""Plot native-Thr similarity distributions for original Tier A/B frames."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LOGGER = logging.getLogger(__name__)
FIELDS = (
    ("og1_shift", "OG1 displacement (Å)", (0, 4.5)),
    ("nterm_shift", "N-terminal N displacement (Å)", (0, 4.5)),
    ("psi_delta", "Absolute ψ deviation (degrees)", (0, 180)),
    ("chi1_delta", "Absolute χ1 deviation (degrees)", (0, 180)),
)
COLORS = {"WT": "#008e9a", "broader_mask07": "#be347d"}


def plot(bundle: Path) -> None:
    """Plot recorded distributions and original cutoffs without reclassifying.

    Args:
        bundle: Existing audit bundle containing the source-bound scalar report.
    """
    report_path = bundle / "scalar_gate_report.json"
    report = json.loads(report_path.read_text())
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for dataset, color in COLORS.items():
        source = Path(report["datasets"][dataset]["source_csv"])
        if (
            hashlib.sha256(source.read_bytes()).hexdigest()
            != report["datasets"][dataset]["source_sha256"]
        ):
            raise ValueError("Source scalar CSV hash changed")
        with source.open() as stream:
            rows = [r for r in csv.DictReader(stream) if r["tier"] in {"A", "B"}]
        group = report["datasets"][dataset]["groups"]["all:tier_ab"]
        if len(rows) != group["n"]:
            raise ValueError("Tier A/B denominator differs from saved audit")
        label = "WT" if dataset == "WT" else "mask_07"
        for ax, (field, title, limits) in zip(axes.flat, FIELDS, strict=True):
            data = np.abs([float(r[field]) for r in rows])
            if not np.isfinite(data).all():
                raise ValueError("Nonfinite input scalar")
            bins = np.linspace(*limits, 51)
            ax.hist(
                data,
                bins=bins,
                density=True,
                histtype="stepfilled",
                alpha=0.16,
                color=color,
            )
            ax.hist(
                data,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=1.5,
                color=color,
                label=f"{label}, n={len(rows):,}",
            )
            ax.set(xlim=limits, xlabel=title, ylabel="Density")
            y = 0.9 if dataset == "WT" else 0.80
            ax.text(
                0.98,
                y,
                f"{label}: {group['fail_percent'][field]:.1f}% fail",
                transform=ax.transAxes,
                ha="right",
                color=color,
                fontsize=10,
            )
    for ax, (field, _, _) in zip(axes.flat, FIELDS, strict=True):
        limit = report["original_limits"][field]
        ax.axvline(limit, color="black", linestyle="--", linewidth=1.3)
        ax.set_title(f"Original cutoff: {limit:g}{' Å' if 'shift' in field else '°'}")
        ax.spines[["top", "right"]].set_visible(False)
    axes.flat[0].legend(loc="upper left", frameon=False, fontsize=9)
    fig.suptitle("Native-Thr similarity among Tier A/B MD frames", fontsize=14)
    fig.text(
        0.5,
        0.015,
        "One replica per sequence; four nested sites, 100-ps sampling. Cutoffs unchanged; not an activity classifier.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    for extension in ("png", "pdf"):
        fig.savefig(bundle / f"gate_distributions.{extension}", dpi=180)
    plt.close(fig)
    LOGGER.info("Saved source-validated native-Thr distributions in %s", bundle)


def compose_overlays(bundle: Path) -> None:
    """Caption shared-camera renders using one common whitespace crop.

    Args:
        bundle: Existing rendered molecule bundle with snapshot metadata.
    """
    from PIL import Image

    manifest = json.loads((bundle / "snapshot_manifest.json").read_text())
    entries = {(s["dataset"], s["role"]): s for s in manifest["snapshots"]}
    roles = ("initial", "typical", "closest")
    images, sources = {}, {}
    bounds = []
    for role in roles:
        path = bundle / f"native_thr_{role}_side_by_side.png"
        sources[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        array = np.asarray(Image.open(path).convert("RGB"))
        midpoint = array.shape[1] // 2
        for index, dataset in enumerate(("WT", "broader_mask07")):
            panel = array[:, index * midpoint : (index + 1) * midpoint]
            images[(dataset, role)] = panel
            ys, xs = np.where(np.any(panel < 235, axis=2))
            if not len(xs):
                raise ValueError("Molecule panel is empty")
            bounds.append((xs.min(), ys.min(), xs.max(), ys.max()))
    xmin = max(0, min(b[0] for b in bounds) - 25)
    ymin = max(0, min(b[1] for b in bounds) - 25)
    xmax = min(midpoint, max(b[2] for b in bounds) + 26)
    ymax = min(array.shape[0], max(b[3] for b in bounds) + 26)

    def figure(selected_roles: tuple[str, ...], filename: str) -> None:
        fig, axes = plt.subplots(
            len(selected_roles),
            2,
            squeeze=False,
            figsize=(10, 5.8 * len(selected_roles)),
        )
        for row, role in enumerate(selected_roles):
            for col, dataset in enumerate(("WT", "broader_mask07")):
                item = entries[(dataset, role)]
                ax = axes[row, col]
                ax.imshow(images[(dataset, role)][ymin:ymax, xmin:xmax])
                label = "WT" if dataset == "WT" else "mask_07"
                when = (
                    "minimized start"
                    if role == "initial"
                    else f"{role}, {item['time_ns']:g} ns"
                )
                ax.set_title(f"{label}: {when} (site {item['site']})", fontsize=12)
                ax.axis("off")
                m = item["metrics"]
                text = (
                    f"OG1 shift {m['og1_shift']:.2f} Å; N shift {m['nterm_shift']:.2f} Å\n"
                    f"Δψ {m['psi_delta']:+.1f}°; Δχ1 {m['chi1_delta']:+.1f}°"
                )
                ax.text(
                    0.5,
                    -0.025,
                    text,
                    transform=ax.transAxes,
                    ha="center",
                    va="top",
                    fontsize=10,
                )
        fig.suptitle(
            "Thr233 aligned to n12c.cif — catalytic Thr excluded from fit",
            fontsize=13,
            y=0.995,
        )
        fig.text(
            0.5,
            0.015,
            "Carbon colors: native reference gray; WT cyan; mask_07 magenta; substrate orange. Same camera and scale.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.07 / len(selected_roles), 1, 0.97))
        fig.savefig(bundle / filename, dpi=180)
        plt.close(fig)

    figure(("typical",), "native_thr_typical_annotated.png")
    figure(roles, "native_thr_gallery.png")
    (bundle / "figure_composition.json").write_text(
        json.dumps(
            {
                "source_image_sha256": sources,
                "common_crop_xyxy": [int(x) for x in (xmin, ymin, xmax, ymax)],
                "geometry_modified": False,
                "scaling": "Same camera and common pixel crop for every panel",
            },
            indent=2,
        )
        + "\n"
    )
    LOGGER.info("Saved captioned molecular overlays with common crop")


def main() -> None:
    """Parse the existing bundle path and render the diagnostic distributions."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--compose-overlays", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.compose_overlays:
        compose_overlays(args.bundle.resolve(strict=True))
    else:
        plot(args.bundle.resolve(strict=True))


if __name__ == "__main__":
    main()
