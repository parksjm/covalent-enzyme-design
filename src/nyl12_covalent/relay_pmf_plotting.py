"""Render standardized relay PMF and overlap diagnostics without refitting."""

from __future__ import annotations

import hashlib
import json
import logging
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import numpy as np

from .relay_pmf_plot_panels import plot_neighbors, plot_overlap, plot_pmf

LOGGER = logging.getLogger(__name__)

PMF_ARRAY_FIELDS = (
    "centers",
    "pmf_kcal_mol",
    "standard_error_kcal_mol",
    "lower_95_kcal_mol",
    "upper_95_kcal_mol",
    "selected_frames_per_bin",
    "weight_ess_per_bin",
    "contributing_windows_per_bin",
)


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file.

    Args:
        path: Existing file to hash.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_inputs(
    report: dict[str, Any], uncertainty: dict[str, Any]
) -> dict[str, Any]:
    """Validate the fixed-grid uncertainty payload against its source report.

    Args:
        report: Immutable relay umbrella analysis report.
        uncertainty: Pointwise PMF uncertainty report from the refitted MBAR
            bootstrap/asymptotic analysis.

    Returns:
        Validated arrays, state labels, horizon, and physical temperature.
    """
    if not isinstance(report.get("pmf"), dict):
        raise ValueError("Analysis report must contain a PMF")
    pmf = report["pmf"]
    centers = np.asarray(uncertainty["centers"], dtype=float)
    report_centers = np.asarray(pmf["centers"], dtype=float)
    if (
        centers.ndim != 1
        or centers.shape != report_centers.shape
        or not np.isfinite(centers).all()
    ):
        raise ValueError("Uncertainty payload must match the frozen report PMF grid")
    arrays: dict[str, np.ndarray] = {"centers": centers}
    for field in PMF_ARRAY_FIELDS[1:]:
        values = uncertainty[field]
        if len(values) != len(centers):
            raise ValueError(
                f"Uncertainty field {field} must match the report PMF grid"
            )
        try:
            arrays[field] = np.asarray(
                [np.nan if value is None else float(value) for value in values],
                dtype=float,
            )
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"Uncertainty field {field} must be numeric or null"
            ) from error
    bootstrap = uncertainty.get("block_bootstrap")
    if bootstrap is not None:
        support = bootstrap["support_fraction"]
        if len(support) != len(centers):
            raise ValueError(
                "Block-bootstrap support_fraction must match the report PMF grid"
            )
        support_values = np.asarray(
            [np.nan if value is None else float(value) for value in support],
            dtype=float,
        )
        finite_support = np.isfinite(support_values)
        if np.any(
            (support_values[finite_support] < 0) | (support_values[finite_support] > 1)
        ):
            raise ValueError("Block-bootstrap support_fraction must be in [0, 1]")
        arrays["block_bootstrap_support_fraction"] = support_values
        for field in ("lower_95_kcal_mol", "upper_95_kcal_mol"):
            values = bootstrap[field]
            if len(values) != len(centers):
                raise ValueError(
                    f"Block-bootstrap field {field} must match the report PMF grid"
                )
            arrays[f"block_bootstrap_{field}"] = np.asarray(
                [np.nan if value is None else float(value) for value in values],
                dtype=float,
            )
    if not np.allclose(centers, report_centers, rtol=0, atol=1e-12):
        raise ValueError("Uncertainty centers differ from the frozen report PMF grid")
    if uncertainty.get("report_sha256") != report.get("_source_sha256"):
        raise ValueError("Uncertainty report_sha256 does not match the analysis report")
    if (
        uncertainty.get("confidence_method")
        != "PyMBAR histogram asymptotic covariance; decorrelated samples"
    ):
        raise ValueError("Unexpected uncertainty confidence method")
    reference_s = float(uncertainty["reference_s"])
    if not np.isfinite(reference_s) or not np.any(
        np.isclose(centers, reference_s, rtol=0, atol=1e-12)
    ):
        raise ValueError("Uncertainty reference_s must match a fixed histogram center")
    if (
        uncertainty.get("sampling_converged") is not False
        or report.get("converged") is not False
    ):
        raise ValueError("These figures require explicit unconverged sampling status")
    horizon = float(report["checkpoint"]["actual_horizon_ps"])
    temperature = float(report["temperature_K"])
    if (
        "temperature_K" in uncertainty
        and float(uncertainty["temperature_K"]) != temperature
    ):
        raise ValueError("Uncertainty and analysis report temperatures differ")
    if (
        not np.isfinite(horizon)
        or horizon <= 0
        or not np.isfinite(temperature)
        or temperature <= 0
    ):
        raise ValueError("Report horizon and temperature must be positive and finite")
    labels = list(report["duration_ps"])
    block_counts = Counter(
        item["window"]
        for item in report.get("blocks", [])
        if item.get("window") in labels
    )
    per_window_block_counts = [block_counts.get(label, 0) for label in labels]
    if per_window_block_counts and min(per_window_block_counts) == 0:
        raise ValueError("Every plotted window needs an accepted block summary")
    overlap = np.asarray(report["overlap"], dtype=float)
    if overlap.shape != (len(labels), len(labels)) or not np.isfinite(overlap).all():
        raise ValueError(
            "Full MBAR overlap matrix must match the ordered window labels"
        )
    gates = report["gates"]["neighbor_histogram_overlap"]
    pairs = np.asarray(gates["pairs"], dtype=int)
    intersections = np.asarray(gates["intersections"], dtype=float)
    if pairs.ndim != 2 or pairs.shape[1] != 2 or len(pairs) != len(intersections):
        raise ValueError("Neighbor overlap pairs and values must have matching shapes")
    if np.any(pairs < 0) or np.any(pairs >= len(labels)):
        raise ValueError("Neighbor comparison pair references an unknown state")
    return dict(
        arrays=arrays,
        reference_s=reference_s,
        horizon_ps=horizon,
        temperature_K=temperature,
        labels=labels,
        blocks_per_window_range=(
            [min(per_window_block_counts), max(per_window_block_counts)]
            if per_window_block_counts
            else None
        ),
        overlap=overlap,
        pairs=pairs,
        intersections=intersections,
        overlap_threshold=float(gates["minimum"]),
    )


def render_diagnostics(
    report_path: Path,
    uncertainty_path: Path,
    output_dir: Path,
    stem: str,
    interval_method: str,
) -> dict[str, Any]:
    """Render standardized figures and provenance from immutable JSON inputs.

    Args:
        report_path: Existing cumulative analysis JSON.
        uncertainty_path: Existing MBAR uncertainty JSON bound to report hash.
        output_dir: Directory for new figures and provenance.
        stem: Descriptive filename stem for all outputs.
        interval_method: Explicit interval source: ``analytic`` or ``block-bootstrap``.

    Returns:
        Provenance mapping containing inputs, hashes, output paths, and status.
    """
    report_path = report_path.resolve(strict=True)
    uncertainty_path = uncertainty_path.resolve(strict=True)
    report_bytes = report_path.read_bytes()
    uncertainty_bytes = uncertainty_path.read_bytes()
    report = json.loads(report_bytes)
    report["_source_sha256"] = hashlib.sha256(report_bytes).hexdigest()
    uncertainty = json.loads(uncertainty_bytes)
    data = validate_inputs(report, uncertainty)
    if interval_method not in {"analytic", "block-bootstrap"}:
        raise ValueError("interval_method must be analytic or block-bootstrap")
    if interval_method == "block-bootstrap" and "block_bootstrap" not in uncertainty:
        raise ValueError("Block-bootstrap intervals were requested but are absent")
    if interval_method == "block-bootstrap":
        if "block_bootstrap_support_fraction" not in data["arrays"]:
            raise ValueError("Block-bootstrap support fractions are missing")
    output_dir = output_dir.resolve()
    destinations = [
        output_dir / f"{stem}_{suffix}.{extension}"
        for suffix in ("pmf", "mbar_overlap", "neighbor_overlap")
        for extension in ("png", "pdf")
    ] + [output_dir / f"{stem}_provenance.json"]
    if any(path.exists() for path in destinations):
        raise FileExistsError(
            "One or more output artifacts already exist; retain prior files"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    with matplotlib.rc_context(
        {
            "font.size": 18,
            "axes.titlesize": 22,
            "axes.labelsize": 21,
            "xtick.labelsize": 17,
            "ytick.labelsize": 17,
        }
    ):
        outputs = []
        outputs.extend(
            plot_pmf(
                data | {"report": report, "uncertainty": uncertainty},
                output_dir,
                stem,
                report["_source_sha256"],
                interval_method,
            )
        )
        outputs.extend(plot_overlap(data, output_dir, stem))
        outputs.extend(plot_neighbors(data, output_dir, stem))
    code_paths = [
        Path(__file__).resolve(),
        Path(__file__).resolve().with_name("relay_pmf_plot_panels.py"),
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "plot_relay_pmf_diagnostics.py",
        Path(__file__).resolve().parents[2]
        / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/plots/plot_relay_pmf_overlap.py",
    ]
    source_hashes = {
        str(path): sha256_file(path)
        for path in [report_path, uncertainty_path, *code_paths]
        if path.is_file()
    }
    provenance = {
        "report_path": str(report_path),
        "uncertainty_path": str(uncertainty_path),
        "source_sha256": source_hashes,
        "uncertainty_upstream_source_sha256": uncertainty.get("source_sha256", {}),
        "output_paths": [str(path) for path in outputs],
        "horizon_ps": data["horizon_ps"],
        "temperature_K": data["temperature_K"],
        "reference_s": data["reference_s"],
        "confidence_method": uncertainty["confidence_method"],
        "interval_method": interval_method,
        "block_bootstrap": uncertainty.get("block_bootstrap"),
        "blocks_per_window_range": data["blocks_per_window_range"],
        "sampling_converged": False,
        "scientific_inputs_or_gates_changed": False,
    }
    provenance_path = output_dir / f"{stem}_provenance.json"
    provenance["output_paths"].append(str(provenance_path))
    provenance_path.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    LOGGER.info(
        "Rendered %d figures; report SHA-256 %s",
        len(outputs),
        source_hashes[str(report_path)],
    )
    return provenance
