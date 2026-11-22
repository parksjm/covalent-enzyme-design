"""Plot raw and relaxed catalytic geometry for the Nyl12 design campaign."""

from __future__ import annotations

import csv
import json
import logging
import math
from pathlib import Path
from typing import Any

import matplotlib
from matplotlib.lines import Line2D

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOGGER = logging.getLogger(__name__)
VARIANTS = ("WT", *(f"variant_{index:02d}" for index in range(1, 9)))
STAGES = ("raw", "relaxed")
SITES = ("AE", "BF", "CG", "DH")
METRICS = ("d1", "d2", "d3", "a_bd", "herschlag_phi")
LABELS = {
    "d1": "Thr OG1–C54 (Å)",
    "d2": "O49–Phe113 N (Å)",
    "d3": "O49–Asn186 ND2 (Å)",
    "a_bd": "Bürgi–Dunitz angle (°)",
    "herschlag_phi": "Herschlag φ (°; signed)",
}
SITE_COLORS = dict(
    zip(SITES, ("#0072B2", "#D55E00", "#009E73", "#CC79A7"), strict=True)
)


def _read_rows(path: Path, required: set[str]) -> list[dict[str, str]]:
    """Read a nonempty CSV with all required columns."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or ())
        if required - fields:
            raise ValueError(f"{path} missing columns: {sorted(required - fields)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} contains no model rows")
    return rows


def _numeric(
    row: dict[str, str], columns: tuple[str, ...], index: int
) -> dict[str, float]:
    """Parse finite metric values without changing sign or substituting blanks."""
    values: dict[str, float] = {}
    for column in columns:
        try:
            value = float(row[column])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Missing or invalid {column} in row {index}") from exc
        if not math.isfinite(value):
            raise ValueError(f"Non-finite {column} in row {index}")
        values[column] = value
    return values


def _check_duplicates(
    rows: list[dict[str, str]], *, site: bool, variants: tuple[str, ...]
) -> None:
    """Reject duplicate model-stage records, including duplicate occupied sites."""
    seen: set[tuple[str, int, str, str | None]] = set()
    for index, row in enumerate(rows, 2):
        variant, stage = row["variant"], row["stage"]
        if variant not in variants or stage not in STAGES:
            raise ValueError(
                f"Unsupported variant/stage in row {index}: {variant}/{stage}"
            )
        try:
            seed = int(row["seed"])
        except ValueError as exc:
            raise ValueError(f"Invalid seed in row {index}") from exc
        site_name = row["site"] if site else None
        if site and site_name not in SITES:
            raise ValueError(f"Invalid site in row {index}: {site_name}")
        key = (variant, seed, stage, site_name)
        if key in seen:
            raise ValueError(f"Duplicate model/stage record: {key}")
        seen.add(key)


def _stage_counts(
    rows: list[dict[str, str]], planned: int
) -> dict[str, dict[str, int]]:
    """Count unique model seeds by stage and expose missing observations."""
    return {
        stage: {
            "planned_models": planned,
            "observed_models": count,
            "missing_models": max(planned - count, 0),
        }
        for stage in STAGES
        for count in [
            len(
                {
                    (row["variant"], int(row["seed"]))
                    for row in rows
                    if row["stage"] == stage
                }
            )
        ]
    }


def _geometry_plot(
    rows: list[dict[str, str]], output: Path, planned: int, variants: tuple[str, ...]
) -> None:
    figure, axes = plt.subplots(2, 5, figsize=(19, 8), sharex=True)
    offsets = dict(zip(SITES, (-0.24, -0.08, 0.08, 0.24), strict=True))
    xpositions = {variant: index for index, variant in enumerate(variants)}
    for row_index, stage in enumerate(STAGES):
        stage_rows = [row for row in rows if row["stage"] == stage]
        n_models = len({(row["variant"], int(row["seed"])) for row in stage_rows})
        missing = max(planned - n_models, 0)
        for column, metric in enumerate(METRICS):
            axis = axes[row_index, column]
            if not stage_rows:
                axis.text(
                    0.5,
                    0.5,
                    f"{missing}/{planned} models missing",
                    transform=axis.transAxes,
                )
            else:
                for row in stage_rows:
                    x = xpositions[row["variant"]] + offsets[row["site"]]
                    value = float(row[metric])
                    if not math.isfinite(value):
                        raise ValueError(f"Non-finite {metric} in geometry plot")
                    axis.scatter(
                        x, value, color=SITE_COLORS[row["site"]], s=22, alpha=0.75
                    )
            axis.set_title(LABELS[metric])
            axis.grid(axis="y", alpha=0.2)
            if column == 0:
                axis.set_ylabel(
                    f"{stage} stage\n{LABELS[metric]}\n{n_models}/{planned} models; {missing} missing"
                )
            if row_index == 1:
                axis.set_xticks(range(len(variants)), variants, rotation=45, ha="right")
    handles = [
        Line2D([], [], marker="o", linestyle="", color=SITE_COLORS[s], label=s)
        for s in SITES
    ]
    figure.legend(
        handles=handles, title="Catalytic site", loc="upper center",
        bbox_to_anchor=(0.5, 0.94), ncol=4,
    )
    figure.suptitle("Occupied-site catalytic geometry by design and stage")
    figure.tight_layout(rect=(0, 0, 1, 0.86))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def _fraction_plot(
    rows: list[dict[str, str]], output: Path, planned: int, variants: tuple[str, ...]
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True)
    xpositions = {variant: index for index, variant in enumerate(variants)}
    stage_style = {"raw": ("#0072B2", -0.12), "relaxed": ("#D55E00", 0.12)}
    for axis, metric, title in zip(
        axes,
        ("tier_ab_fraction", "joint_pass_fraction"),
        ("Model-level Tier A/B fraction", "Model-level joint-pass fraction"),
        strict=True,
    ):
        metric_rows = rows
        counts = _stage_counts(metric_rows, planned)
        raw = counts["raw"]["observed_models"]
        relaxed = counts["relaxed"]["observed_models"]
        raw_missing = counts["raw"]["missing_models"]
        relaxed_missing = counts["relaxed"]["missing_models"]
        axis.set_title(
            f"{title}\nraw {raw}/{planned} ({raw_missing} missing); "
            f"relaxed {relaxed}/{planned} ({relaxed_missing} missing)"
        )
        axis.set_ylim(-0.04, 1.04)
        axis.set_ylabel("Fraction of occupied sites")
        axis.grid(axis="y", alpha=0.2)
        for row in metric_rows:
            value = _numeric(row, (metric,), 0)[metric]
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{metric} must lie in [0, 1]")
            color, offset = stage_style[row["stage"]]
            axis.scatter(
                xpositions[row["variant"]] + offset,
                value,
                color=color,
                s=30,
                alpha=0.75,
            )
        axis.set_xticks(range(len(variants)), variants, rotation=45, ha="right")
    handles = [
        Line2D([], [], marker="o", linestyle="", color=color, label=stage)
        for stage, (color, _) in stage_style.items()
    ]
    figure.legend(
        handles=handles, title="Stage", loc="upper center",
        bbox_to_anchor=(0.5, 0.92), ncol=2,
    )
    figure.suptitle("Model-level geometry pass fractions (sites nested within models)")
    figure.tight_layout(rect=(0, 0, 1, 0.82))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_comparison(
    occupied_csv: Path,
    model_summary_csv: Path,
    output: Path,
    planned_models: int = 45,
    variant_order: tuple[str, ...] = VARIANTS,
) -> dict[str, Any]:
    """Write geometry and model-level pass-fraction plots.

    Args:
        occupied_csv: Per-model, per-site occupied geometry table.
        model_summary_csv: Model-level tier and joint-pass fractions.
        output: Existing or new directory for PNGs and plot manifest.
        planned_models: Expected unique model seeds per stage.
        variant_order: Explicit allowed cohort IDs in plotting order.
    Returns:
        Plot paths and observed/missing model counts.
    """
    if planned_models < 1:
        raise ValueError("planned_models must be positive")
    if (
        not variant_order
        or any(not isinstance(name, str) or not name.strip() for name in variant_order)
        or len(set(variant_order)) != len(variant_order)
    ):
        raise ValueError("variant_order must contain unique nonempty cohort IDs")
    occupied = _read_rows(occupied_csv, {"variant", "seed", "stage", "site", *METRICS})
    summary = _read_rows(
        model_summary_csv,
        {
            "variant",
            "seed",
            "stage",
            "tier_ab_fraction",
            "joint_pass_fraction",
        },
    )
    _check_duplicates(occupied, site=True, variants=variant_order)
    _check_duplicates(summary, site=False, variants=variant_order)
    output.mkdir(parents=True, exist_ok=True)
    geometry_path = output / "geometry_before_after.png"
    fraction_path = output / "tier_native_comparison.png"
    _geometry_plot(occupied, geometry_path, planned_models, variant_order)
    _fraction_plot(summary, fraction_path, planned_models, variant_order)
    manifest: dict[str, Any] = {
        "geometry_plot": str(geometry_path.resolve()),
        "fraction_plot": str(fraction_path.resolve()),
        "planned_models_per_stage": planned_models,
        "variant_order": list(variant_order),
        "occupied_geometry_counts": _stage_counts(occupied, planned_models),
        "model_fraction_counts": _stage_counts(summary, planned_models),
        "replicate_unit": "model seed; sites within a model are nested, not independent",
        "signed_herschlag_phi": True,
        "activity_ranking_or_statistical_test": False,
    }
    (output / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Wrote design comparison plots to %s", output)
    return manifest
