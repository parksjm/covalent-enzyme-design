"""Tests for matched relay PMF comparison plotting."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent.relay_pmf_comparison import (
    compare_reports,
    validate_matched_reports,
)


def _payload(policy: str | None) -> dict[str, object]:
    centers = [0.0, 0.25, 0.5, 0.75]
    payload: dict[str, object] = {
        "centers": centers,
        "pmf_kcal_mol": [0.0, 1.0, 2.0, 1.0],
        "reference_s": 0.0,
        "temperature_K": 300.0,
        "report_sha256": "a" * 64,
        "request_sha256": "c" * 64,
        "lower_95_kcal_mol": [-0.2, 0.6, 1.4, 0.5],
        "upper_95_kcal_mol": [0.2, 1.4, 2.6, 1.5],
        "sampling_converged": False,
        "block_bootstrap": {
            "lower_95_kcal_mol": [-0.2, 0.6, 1.4, 0.5],
            "upper_95_kcal_mol": [0.2, 1.4, 2.6, 1.5],
            "support_fraction": [1.0, 1.0, 1.0, 0.0],
            "reference_s": 0.0,
            "replicates": 200,
        },
    }
    if policy is not None:
        payload["selection_policy"] = policy
    if policy == "continuous_window":
        payload["window_selection"] = [dict(window="w0", duration_ps=20)]
        payload["parent_report_sha256"] = payload["report_sha256"]
        payload["report_sha256"] = "d" * 64
        payload["window_bootstrap"] = payload.pop("block_bootstrap")
    return payload


def test_validate_matched_reports_stops_on_grid_reference_and_source_mismatch() -> None:
    baseline = _payload("blockwise_window")
    continuous = _payload("continuous_window")
    result = validate_matched_reports(baseline, continuous, "analytic")
    assert result["centers"].shape == (4,)
    assert result["temperature_K"] == 300.0
    assert result["horizon_ps"] == 20.0
    for field, value, message in (
        ("centers", [0, 0.2, 0.5, 0.75], "histogram centers"),
        ("reference_s", 0.25, "reference_s"),
        ("parent_report_sha256", "b" * 64, "parent_report_sha256"),
    ):
        altered = _payload("continuous_window")
        altered[field] = value
        with pytest.raises(ValueError, match=message):
            validate_matched_reports(baseline, altered, "analytic")


def test_compare_reports_saves_png_pdf_and_provenance(tmp_path: Path) -> None:
    baseline_path = tmp_path / "baseline.json"
    continuous_path = tmp_path / "continuous.json"
    baseline_path.write_text(json.dumps(_payload(None)))
    continuous_path.write_text(json.dumps(_payload("continuous_window")))
    stem = tmp_path / "comparison_20ps_continuous_vs_blockwise_20261126"
    provenance = compare_reports(baseline_path, continuous_path, stem, "analytic")
    assert all(Path(path).is_file() for path in provenance["output_paths"])
    assert provenance["report_sha256"] == "a" * 64
    assert provenance["source_sha256"][str(baseline_path.resolve())]
    assert provenance["source_sha256"][str(continuous_path.resolve())]
    assert "not independent replicates" in provenance["interpretation"]
    with pytest.raises(FileExistsError, match="already exist"):
        compare_reports(baseline_path, continuous_path, stem, "analytic")
