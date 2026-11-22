"""Tests for incomplete-campaign geometry and native-tier plots."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from nyl12_covalent.design_geometry_plots import _numeric, plot_comparison

GEOMETRY_FIELDS = [
    "variant",
    "seed",
    "stage",
    "site",
    "d1",
    "d2",
    "d3",
    "a_bd",
    "herschlag_phi",
]
SUMMARY_FIELDS = ["variant", "seed", "stage", "tier_ab_fraction", "joint_pass_fraction"]


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    geometry = tmp_path / "occupied.csv"
    summary = tmp_path / "summary.csv"
    geometry_rows = [
        {
            "variant": "WT",
            "seed": "2026112201",
            "stage": "raw",
            "site": site,
            "d1": "2.8",
            "d2": "2.9",
            "d3": "2.7",
            "a_bd": "94.0",
            "herschlag_phi": "-120.0",
        }
        for site in ("AE", "BF", "CG", "DH")
    ]
    _write_csv(geometry, GEOMETRY_FIELDS, geometry_rows)
    _write_csv(
        summary,
        SUMMARY_FIELDS,
        [
            {
                "variant": "WT",
                "seed": "2026112201",
                "stage": "raw",
                "tier_ab_fraction": "0.5",
                "joint_pass_fraction": "0.25",
            }
        ],
    )
    return geometry, summary


def test_partial_campaign_outputs_plots_with_missing_counts_and_signed_phi(
    tmp_path: Path,
) -> None:
    geometry, summary = _inputs(tmp_path)
    output = tmp_path / "plots"
    manifest = plot_comparison(geometry, summary, output, planned_models=45)
    assert manifest["occupied_geometry_counts"] == {
        "raw": {"planned_models": 45, "observed_models": 1, "missing_models": 44},
        "relaxed": {"planned_models": 45, "observed_models": 0, "missing_models": 45},
    }
    assert manifest["signed_herschlag_phi"] is True
    assert manifest["model_fraction_counts"]["relaxed"]["missing_models"] == 45
    for name in ("geometry_before_after.png", "tier_native_comparison.png"):
        path = output / name
        assert path.is_file() and path.stat().st_size > 1000
        assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert (
        json.loads((output / "plot_manifest.json").read_text())[
            "activity_ranking_or_statistical_test"
        ]
        is False
    )


@pytest.mark.parametrize("problem", ["missing_metric", "duplicate", "nan"])
def test_invalid_geometry_rows_fail_explicitly(tmp_path: Path, problem: str) -> None:
    geometry, summary = _inputs(tmp_path)
    rows = list(csv.DictReader(geometry.open(newline="")))
    if problem == "missing_metric":
        fields = [field for field in GEOMETRY_FIELDS if field != "d2"]
        _write_csv(
            geometry,
            fields,
            [{key: value for key, value in rows[0].items() if key != "d2"}],
        )
    elif problem == "duplicate":
        _write_csv(geometry, GEOMETRY_FIELDS, [rows[0], rows[0]])
    else:
        rows[0]["d1"] = "nan"
        _write_csv(geometry, GEOMETRY_FIELDS, rows[:1])
    with pytest.raises(ValueError):
        plot_comparison(geometry, summary, tmp_path / "out")


def test_invalid_fraction_is_not_clamped(tmp_path: Path) -> None:
    geometry, summary = _inputs(tmp_path)
    rows = list(csv.DictReader(summary.open(newline="")))
    rows[0]["joint_pass_fraction"] = "1.1"
    _write_csv(summary, SUMMARY_FIELDS, rows)
    with pytest.raises(ValueError, match=r"must lie in \[0, 1\]"):
        plot_comparison(geometry, summary, tmp_path / "out")


def test_numeric_parser_preserves_signed_phi() -> None:
    assert _numeric({"herschlag_phi": "-120"}, ("herschlag_phi",), 2) == {
        "herschlag_phi": -120.0
    }


def test_explicit_mask_and_control_cohort_ids_are_supported(tmp_path: Path) -> None:
    geometry, summary = _inputs(tmp_path)
    geometry_rows = list(csv.DictReader(geometry.open(newline="")))
    summary_rows = list(csv.DictReader(summary.open(newline="")))
    variants = ("WT", "control_K59D", "mask_01")
    _write_csv(
        geometry,
        GEOMETRY_FIELDS,
        [{**row, "variant": name} for name in variants for row in geometry_rows],
    )
    _write_csv(
        summary,
        SUMMARY_FIELDS,
        [{**row, "variant": name} for name in variants for row in summary_rows],
    )
    manifest = plot_comparison(
        geometry, summary, tmp_path / "plots", planned_models=3, variant_order=variants
    )
    assert manifest["variant_order"] == list(variants)
    assert manifest["model_fraction_counts"]["raw"]["observed_models"] == 3
    with pytest.raises(ValueError, match="Unsupported variant/stage"):
        plot_comparison(geometry, summary, tmp_path / "undeclared")


@pytest.mark.parametrize("variants", [(), ("WT", "WT"), ("",), (" ",)])
def test_invalid_cohort_order_fails_before_plotting(
    tmp_path: Path, variants: tuple[str, ...]
) -> None:
    geometry, summary = _inputs(tmp_path)
    with pytest.raises(ValueError, match="unique nonempty cohort IDs"):
        plot_comparison(geometry, summary, tmp_path / "plots", variant_order=variants)
