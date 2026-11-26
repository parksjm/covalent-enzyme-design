"""Tests for standardized report-only relay PMF diagnostics."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.relay_pmf_plotting import render_diagnostics, validate_inputs


def _inputs(tmp_path: Path) -> tuple[Path, Path, dict[str, object]]:
    centers = np.linspace(-0.1, 1.1, 60).tolist()
    labels = [f"state_{index:02d}" for index in range(57)]
    report: dict[str, object] = {
        "pmf": {"centers": centers, "edges": np.linspace(-0.11, 1.11, 61).tolist()},
        "checkpoint": {"actual_horizon_ps": 20.0},
        "temperature_K": 300.0,
        "converged": False,
        "duration_ps": dict.fromkeys(labels, 20.0),
        "blocks": [
            {"window": label, "raw_frames": 100} for label in labels for _ in range(4)
        ],
        "overlap": (
            np.eye(57) * 0.9 + (np.ones((57, 57)) - np.eye(57)) * 0.1 / 56
        ).tolist(),
        "gates": {
            "neighbor_histogram_overlap": {
                "pairs": [[index, index + 1] for index in range(56)] + [[29, 31]],
                "intersections": [0.25] * 56 + [0.41],
                "minimum": 0.2,
            }
        },
    }
    report_path = tmp_path / "analysis.json"
    report_path.write_text(json.dumps(report))
    report_sha = hashlib.sha256(report_path.read_bytes()).hexdigest()
    values: list[float | None] = [0.0] * 60
    values[21] = None
    uncertainty = {
        "centers": centers,
        "pmf_kcal_mol": values,
        "standard_error_kcal_mol": [0.1] * 60,
        "lower_95_kcal_mol": [-0.2] * 60,
        "upper_95_kcal_mol": [0.2] * 60,
        "selected_frames_per_bin": [10] * 60,
        "weight_ess_per_bin": [8.0] * 60,
        "contributing_windows_per_bin": [3] * 60,
        "reference_s": centers[10],
        "confidence_method": "PyMBAR histogram asymptotic covariance; decorrelated samples",
        "block_bootstrap": {
            "lower_95_kcal_mol": [-0.3] * 60,
            "upper_95_kcal_mol": [0.3] * 60,
            "support_fraction": [1.0] * 60,
            "replicates": 200,
            "method": "within-window 5-ps moving-block bootstrap; MBAR refit",
        },
        "sampling_converged": False,
        "report_sha256": report_sha,
    }
    uncertainty_path = tmp_path / "uncertainty.json"
    uncertainty_path.write_text(json.dumps(uncertainty))
    return report_path, uncertainty_path, report


def test_validate_inputs_requires_fixed_grid_and_matching_report_hash(
    tmp_path: Path,
) -> None:
    report_path, uncertainty_path, report = _inputs(tmp_path)
    report["_source_sha256"] = hashlib.sha256(report_path.read_bytes()).hexdigest()
    uncertainty = json.loads(uncertainty_path.read_text())
    validated = validate_inputs(report, uncertainty)
    assert len(validated["arrays"]["centers"]) == 60
    assert validated["labels"][-1] == "state_56"
    uncertainty["report_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="report_sha256"):
        validate_inputs(report, uncertainty)


def test_render_diagnostics_writes_png_pdf_and_hash_provenance(tmp_path: Path) -> None:
    report_path, uncertainty_path, _ = _inputs(tmp_path)
    output_dir = tmp_path / "figures"
    provenance = render_diagnostics(
        report_path,
        uncertainty_path,
        output_dir,
        "relay_20ps",
        "block-bootstrap",
    )
    assert len(provenance["output_paths"]) == 7
    assert all(Path(path).is_file() for path in provenance["output_paths"])
    record = json.loads((output_dir / "relay_20ps_provenance.json").read_text())
    assert record["horizon_ps"] == 20.0
    assert record["temperature_K"] == 300.0
    assert record["sampling_converged"] is False
    assert record["interval_method"] == "block-bootstrap"
    assert str(report_path.resolve()) in record["source_sha256"]
    with pytest.raises(FileExistsError, match="already exist"):
        render_diagnostics(
            report_path,
            uncertainty_path,
            output_dir,
            "relay_20ps",
            "block-bootstrap",
        )
