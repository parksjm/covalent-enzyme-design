"""Compare matched relay PMFs from alternate correlation treatments."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LOGGER = logging.getLogger(__name__)


def _digest(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _numeric(values: list[Any]) -> np.ndarray:
    """Convert nullable JSON numeric arrays to NumPy arrays."""
    return np.asarray([np.nan if value is None else float(value) for value in values])


def _arrays(payload: dict[str, Any], interval_method: str) -> dict[str, np.ndarray]:
    """Read pointwise PMF and selected confidence interval arrays."""
    arrays = {
        "centers": _numeric(payload["centers"]),
        "pmf_kcal_mol": _numeric(payload["pmf_kcal_mol"]),
    }
    if interval_method == "analytic":
        arrays["lower_95_kcal_mol"] = _numeric(payload["lower_95_kcal_mol"])
        arrays["upper_95_kcal_mol"] = _numeric(payload["upper_95_kcal_mol"])
        arrays["support_fraction"] = np.isfinite(arrays["pmf_kcal_mol"]).astype(float)
    else:
        policy = (
            "block_bootstrap" if "block_bootstrap" in payload else "window_bootstrap"
        )
        bootstrap = payload[policy]
        arrays["lower_95_kcal_mol"] = _numeric(bootstrap["lower_95_kcal_mol"])
        arrays["upper_95_kcal_mol"] = _numeric(bootstrap["upper_95_kcal_mol"])
        arrays["support_fraction"] = _numeric(bootstrap["support_fraction"])
    return arrays


def validate_matched_reports(
    baseline: dict[str, Any], continuous: dict[str, Any], interval_method: str
) -> dict[str, Any]:
    """Validate matched grids, reference, source, temperature, and intervals.

    Args:
        baseline: Archived blockwise uncertainty report.
        continuous: New continuous-window uncertainty report.
        interval_method: Explicit interval source, ``analytic`` or ``bootstrap``.

    Returns:
        Validated arrays and shared comparison metadata.
    """
    if continuous.get("selection_policy") != "continuous_window":
        raise ValueError(
            "Continuous report must declare selection_policy=continuous_window"
        )
    durations = [float(row["duration_ps"]) for row in continuous["window_selection"]]
    if (
        not durations
        or not np.isfinite(durations).all()
        or min(durations) <= 0
        or not np.allclose(durations, durations[0])
    ):
        raise ValueError("Continuous windows must share one positive sampled horizon")
    if interval_method not in {"analytic", "bootstrap"}:
        raise ValueError("interval_method must be analytic or bootstrap")
    if continuous.get("parent_report_sha256") != baseline.get("report_sha256"):
        raise ValueError(
            "Continuous parent_report_sha256 must identify the baseline report"
        )
    for field in ("request_sha256", "temperature_K", "reference_s"):
        if baseline.get(field) != continuous.get(field):
            raise ValueError(f"Cannot compare reports with different {field}")
    for payload in (baseline, continuous):
        bootstrap = payload.get("block_bootstrap", payload.get("window_bootstrap", {}))
        nested_reference = bootstrap.get("reference_s")
        if nested_reference is not None and not np.isclose(
            float(nested_reference), float(payload["reference_s"]), rtol=0, atol=1e-12
        ):
            raise ValueError("Interval does not use the declared common reference")
    baseline_arrays = _arrays(baseline, interval_method)
    continuous_arrays = _arrays(continuous, interval_method)
    for label, arrays, payload in (
        ("baseline", baseline_arrays, baseline),
        ("continuous", continuous_arrays, continuous),
    ):
        centers = arrays["centers"]
        if centers.ndim != 1 or not len(centers) or not np.isfinite(centers).all():
            raise ValueError(f"{label} centers must be finite and one-dimensional")
        if any(values.shape != centers.shape for values in arrays.values()):
            raise ValueError(f"{label} uncertainty arrays do not match centers")
        support = arrays["support_fraction"]
        if np.any(np.isfinite(support) & ((support < 0) | (support > 1))):
            raise ValueError(f"{label} support_fraction must be in [0, 1]")
        if payload.get("sampling_converged") is not False:
            raise ValueError(f"{label} report must retain explicit unconverged status")
    if not np.allclose(
        baseline_arrays["centers"], continuous_arrays["centers"], rtol=0, atol=1e-12
    ):
        raise ValueError("Cannot compare reports with different histogram centers")
    return dict(
        centers=baseline_arrays["centers"],
        baseline=baseline_arrays,
        continuous=continuous_arrays,
        temperature_K=float(baseline["temperature_K"]),
        reference_s=float(baseline["reference_s"]),
        horizon_ps=durations[0],
    )


def compare_reports(
    baseline_path: Path,
    continuous_path: Path,
    output_stem: Path,
    interval_method: str,
) -> dict[str, Any]:
    """Render matched PMFs and retain source hashes in provenance.

    Args:
        baseline_path: Archived blockwise uncertainty JSON.
        continuous_path: Continuous-window uncertainty JSON.
        output_stem: Full output path stem, without a file extension.
        interval_method: Explicit interval source, ``analytic`` or ``bootstrap``.

    Returns:
        Provenance metadata for the PNG, PDF, and JSON outputs.
    """
    baseline_path = baseline_path.resolve(strict=True)
    continuous_path = continuous_path.resolve(strict=True)
    baseline_bytes = baseline_path.read_bytes()
    continuous_bytes = continuous_path.read_bytes()
    baseline = json.loads(baseline_bytes)
    continuous = json.loads(continuous_bytes)
    data = validate_matched_reports(baseline, continuous, interval_method)
    output_stem = output_stem.resolve()
    destinations = [
        output_stem.with_suffix(suffix) for suffix in (".png", ".pdf", ".json")
    ]
    if any(path.exists() for path in destinations):
        raise FileExistsError("Comparison outputs already exist; retain prior files")
    output_stem.parent.mkdir(parents=True, exist_ok=True)

    colors = {"blockwise": "#2366a8", "continuous": "#e67e22", "reference": "#707070"}
    centers = data["centers"]
    figure, (axis, delta_axis) = plt.subplots(
        2,
        1,
        figsize=(14, 12),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1]},
        facecolor="white",
    )
    for name, color, label in (
        ("baseline", colors["blockwise"], "Blockwise 5-ps selection"),
        ("continuous", colors["continuous"], "Continuous-window selection"),
    ):
        arrays = data[name]
        pmf = arrays["pmf_kcal_mol"]
        lower = arrays["lower_95_kcal_mol"]
        upper = arrays["upper_95_kcal_mol"]
        support = arrays["support_fraction"]
        valid = np.isfinite(pmf)
        interval = valid & np.isfinite(lower) & np.isfinite(upper) & (support >= 1)
        axis.plot(
            centers,
            np.where(valid, pmf, np.nan),
            "o-",
            lw=2.2,
            ms=5,
            color=color,
            label=label,
        )
        interval_label = (
            "matched asymptotic 95% CI"
            if interval_method == "analytic"
            else "pointwise 95% bootstrap CI"
        )
        axis.fill_between(
            centers,
            lower,
            upper,
            where=interval,
            interpolate=False,
            color=color,
            alpha=0.20,
            label=f"{label}: {interval_label}",
        )
    axis.plot(
        data["reference_s"],
        0,
        "D",
        color=colors["reference"],
        ms=10,
        label="Common E–S reference",
    )
    axis.axhline(0, color=colors["reference"], lw=1, ls="--")
    axis.set_title(
        f"Matched relay PMFs · {data['temperature_K']:g} K · {data['horizon_ps']:g} ps/window · sampling unconverged",
        fontsize=22,
    )
    axis.set_ylabel("PMF relative to common E–S reference\n(kcal mol⁻¹)", fontsize=21)
    caption = (
        "Matched asymptotic bands; statistical uncertainty only"
        if interval_method == "analytic"
        else "Bootstrap units differ: archived5psblocks versus whole-window-decorrelated configurations"
    )
    axis.text(
        0.015,
        0.97,
        "Shared trajectory data; intervals are not independent-sample significance tests\n"
        + caption,
        transform=axis.transAxes,
        va="top",
        fontsize=14,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.92},
    )
    axis.legend(loc="best", fontsize=16)
    axis.tick_params(labelsize=17)
    delta_valid = np.isfinite(data["baseline"]["pmf_kcal_mol"]) & np.isfinite(
        data["continuous"]["pmf_kcal_mol"]
    )
    delta = data["continuous"]["pmf_kcal_mol"] - data["baseline"]["pmf_kcal_mol"]
    delta_axis.plot(
        centers, np.where(delta_valid, delta, np.nan), "o-", color="#444444", lw=1.8
    )
    delta_axis.axhline(0, color=colors["reference"], lw=1, ls="--")
    delta_axis.set_ylabel("Continuous − blockwise\n(kcal mol⁻¹)", fontsize=21)
    delta_axis.set_xlabel("Normalized path coordinate, $s$", fontsize=21)
    delta_axis.tick_params(labelsize=17)
    for panel in (axis, delta_axis):
        panel.grid(axis="y", alpha=0.2)
        panel.set_axisbelow(True)
    figure.savefig(destinations[0], dpi=200)
    figure.savefig(destinations[1])
    plt.close(figure)

    source_hashes = {
        str(baseline_path): hashlib.sha256(baseline_bytes).hexdigest(),
        str(continuous_path): hashlib.sha256(continuous_bytes).hexdigest(),
    }
    for source in (
        Path(__file__).resolve(),
        Path(__file__).resolve().parents[2] / "scripts/compare_relay_pmf_selection.py",
    ):
        if source.is_file():
            source_hashes[str(source)] = _digest(source)
    provenance = {
        "baseline_path": str(baseline_path),
        "continuous_path": str(continuous_path),
        "baseline_selection_policy": baseline.get(
            "selection_policy", "archived blockwise window selection"
        ),
        "continuous_selection_policy": continuous["selection_policy"],
        "report_sha256": baseline["report_sha256"],
        "request_sha256": baseline["request_sha256"],
        "reference_s": data["reference_s"],
        "temperature_K": data["temperature_K"],
        "horizon_ps": data["horizon_ps"],
        "interval_method": interval_method,
        "source_sha256": source_hashes,
        "baseline_upstream_source_sha256": baseline.get("source_sha256", {}),
        "continuous_upstream_source_sha256": continuous.get("source_sha256", {}),
        "output_paths": [str(path) for path in destinations],
        "interpretation": "Paired descriptive comparison of shared trajectory data; not independent replicates or a significance test.",
    }
    destinations[2].write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    LOGGER.info("Saved matched PMF comparison to %s", output_stem)
    return provenance
