"""Plots for retrospective filter calibration reports."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_LOG = logging.getLogger(__name__)
_GATES = ("d1", "d2", "d3", "a1", "abs_phi")


def _validated(report: dict[str, Any]) -> list[dict[str, Any]]:
    cohorts = report["cohorts"]
    if not cohorts:
        raise ValueError("Calibration report must contain at least one cohort")
    for cohort in cohorts:
        for variant in cohort["per_variant"]:
            n = variant["n_models"]
            if not isinstance(n, int) or n <= 0:
                raise ValueError(
                    f"{cohort['id']}/{variant['variant']}: n_models must be positive"
                )
            gates = variant["gate_fail_counts"]
            missing = set(_GATES) - gates.keys()
            if missing:
                raise ValueError(f"Missing required gate counts: {sorted(missing)}")
            counts = [gates[key] for key in _GATES]
            if any(
                not isinstance(value, int) or not 0 <= value <= n for value in counts
            ):
                raise ValueError(
                    f"{cohort['id']}/{variant['variant']}: invalid gate count"
                )
            for key in ("n_tier_a", "n_tier_ab"):
                if not isinstance(variant[key], int) or not 0 <= variant[key] <= n:
                    raise ValueError(
                        f"{cohort['id']}/{variant['variant']}: invalid {key}"
                    )
            if variant["n_tier_a"] > variant["n_tier_ab"]:
                raise ValueError(
                    f"{cohort['id']}/{variant['variant']}: n_tier_a exceeds n_tier_ab"
                )
    return cohorts


def _calibration_figure(cohorts: list[dict[str, Any]]) -> plt.Figure:
    fig, axes = plt.subplots(
        len(cohorts), 2, figsize=(12, 4.2 * len(cohorts)), squeeze=False
    )
    for row, cohort in enumerate(cohorts):
        variants = cohort["per_variant"]
        labels = [item["variant"] for item in variants]
        x = np.arange(len(variants))
        tier_a = np.array([item["n_tier_a"] for item in variants], dtype=float)
        tier_ab = np.array([item["n_tier_ab"] for item in variants], dtype=float)
        n_models = np.array([item["n_models"] for item in variants], dtype=float)
        exclusive_b = tier_ab - tier_a
        ax = axes[row, 0]
        ax.bar(x, tier_a / n_models, label="Tier A")
        ax.bar(
            x,
            exclusive_b / n_models,
            bottom=tier_a / n_models,
            label="Exclusive Tier B",
        )
        ax.set_xticks(x, labels, rotation=30, ha="right")
        ax.set_ylim(0, 1.18)
        ax.set_ylabel("Fraction of models")
        ax.set_title(
            f"{cohort['id']} · {cohort['system']} · {cohort['register']}\nTier retention"
        )
        ax.legend(frameon=False)
        for idx, item in enumerate(variants):
            ax.text(
                idx,
                min(tier_ab[idx] / n_models[idx] + 0.025, 1.04),
                f"{int(tier_ab[idx])}/{int(n_models[idx])}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

        heat = np.array(
            [
                [item["gate_fail_counts"][gate] / item["n_models"] for gate in _GATES]
                for item in variants
            ],
            dtype=float,
        )
        image = axes[row, 1].imshow(heat, vmin=0, vmax=1, cmap="viridis", aspect="auto")
        axes[row, 1].set_xticks(np.arange(len(_GATES)), _GATES)
        axes[row, 1].set_yticks(x, labels)
        axes[row, 1].set_title("Per-gate failure fraction")
        fig.colorbar(image, ax=axes[row, 1], label="Failure fraction")
    fig.suptitle(
        "Retrospective static-model diagnostics; no validated activity predictor", y=1.0
    )
    fig.tight_layout()
    return fig


def _activity_figure(cohorts: list[dict[str, Any]]) -> plt.Figure:
    fig, axes = plt.subplots(
        1, len(cohorts), figsize=(6 * len(cohorts), 5.5), squeeze=False
    )
    for idx, cohort in enumerate(cohorts):
        ax = axes[0, idx]
        points = [
            item
            for item in cohort["per_variant"]
            if item["selected_endpoint_relative_activity"] is not None
        ]
        x = [item["tier_ab_fraction"] for item in points]
        y = [item["selected_endpoint_relative_activity"] for item in points]
        for item, px, py in zip(points, x, y):
            ax.scatter(
                px, py, marker="s" if item["variant"].upper() == "WT" else "o", s=55
            )
        for item, px, py in zip(points, x, y):
            ax.annotate(
                item["variant"],
                (px, py),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
            )
            if py == 0:
                ax.annotate(
                    "zero product signal",
                    (px, py),
                    xytext=(4, -12),
                    textcoords="offset points",
                    fontsize=7,
                )
        ax.set_xlabel("Tier AB fraction")
        ax.set_ylabel("Reported product signal (WT-relative fold)")
        ax.set_title(f"{cohort['id']} · {cohort['system']} · {cohort['register']}")
        ax.scatter([], [], marker="s", color="black", label="WT")
        ax.scatter([], [], marker="o", color="black", label="Other variants")
        ax.legend(frameon=False)
    fig.suptitle("Polymer-assay product signal vs modeled-register retention")
    fig.tight_layout()
    return fig


def export_calibration(report: dict[str, Any], output: Path) -> None:
    """Export calibration and activity-versus-retention diagnostic plots.

    Args:
        report: Calibration report containing separate cohort summaries.
        output: Directory receiving ``calibration.png``, ``calibration.svg``,
            and ``activity_vs_retention.png``.

    Raises:
        KeyError: If required report fields are missing.
        ValueError: If required gate counts or model denominators are invalid.
    """
    cohorts = _validated(report)
    output.mkdir(parents=True, exist_ok=True)
    fig = _calibration_figure(cohorts)
    fig.savefig(output / "calibration.png", dpi=180, bbox_inches="tight")
    fig.savefig(output / "calibration.svg", bbox_inches="tight")
    plt.close(fig)
    fig = _activity_figure(cohorts)
    fig.savefig(output / "activity_vs_retention.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    _LOG.info("Wrote calibration plots to %s", output)
