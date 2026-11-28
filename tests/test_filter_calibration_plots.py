"""Smoke tests for calibration plot export."""

from pathlib import Path

import pytest

from nyl12_covalent.filter_calibration_plots import export_calibration


def _report() -> dict:
    variants = []
    for name, tier_a, tier_ab, signal in (("WT", 2, 3, 1.0), ("L4", 1, 2, 0.0)):
        variants.append(
            {
                "variant": name,
                "n_models": 4,
                "n_tier_a": tier_a,
                "n_tier_ab": tier_ab,
                "tier_ab_fraction": tier_ab / 4,
                "gate_fail_counts": {"d1": 1, "d2": 2, "d3": 0, "a1": 1, "abs_phi": 2},
                "gate_only_failure_counts": {},
                "n_phi_positive": 1,
                "n_phi_negative": 1,
                "n_phi_zero": 2,
                "selected_endpoint_relative_activity": signal,
                "endpoint_label": "product",
                "any_product_reported_positive": signal > 0,
                "median_interface_delta_tier_ab": None,
            }
        )
    return {
        "cohorts": [
            {
                "id": "PA6",
                "system": "PA6",
                "register": "r1",
                "experimental_comparison": {},
                "per_variant": variants,
                "summary": {},
                "per_model": [],
            }
        ]
    }


def test_export_writes_named_nonempty_figures(tmp_path: Path) -> None:
    export_calibration(_report(), tmp_path)
    for name in ("calibration.png", "calibration.svg", "activity_vs_retention.png"):
        assert (tmp_path / name).stat().st_size > 0


def test_export_rejects_missing_gate(tmp_path: Path) -> None:
    report = _report()
    del report["cohorts"][0]["per_variant"][0]["gate_fail_counts"]["abs_phi"]
    with pytest.raises(ValueError, match="Missing required gate counts"):
        export_calibration(report, tmp_path)


def test_export_rejects_tier_a_above_tier_ab(tmp_path: Path) -> None:
    report = _report()
    report["cohorts"][0]["per_variant"][0]["n_tier_a"] = 4
    with pytest.raises(ValueError, match="n_tier_a exceeds n_tier_ab"):
        export_calibration(report, tmp_path)
