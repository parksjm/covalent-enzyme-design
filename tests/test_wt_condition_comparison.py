"""Tests for explicit positive-phi WT condition comparisons."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from design.nyl12_relay_ts_ensemble_20261126.broader_pilot.wt_paired_precursor_controls_20261127.campaign.wt_condition_comparison import (
    build_condition_comparison,
)


def _rows() -> list[dict]:
    rows = []
    for site in ("AE", "BF", "CG", "DH"):
        for ligand in ("I", "J", "K", "L"):
            for amide in range(3):
                rows.append(
                    {
                        "site": site,
                        "ligand": ligand,
                        "carbonyl_c": f"C{amide}",
                        "carbonyl_o": f"O{amide}",
                        "leaving_n": f"N{amide}",
                        "d1": 3.0 + amide,
                        "herschlag_phi": 70.0,
                        "tier": "A" if amide == 0 else "reject",
                    }
                )
    return rows


def _write_metric(path: Path, *, phase: str = "native") -> None:
    rows = _rows()
    assigned = [
        next(
            r
            for r in rows
            if r["site"] == site and r["ligand"] == "I" and r["carbonyl_c"] == "C0"
        )
        for site in ("AE", "BF", "CG", "DH")
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "all_amide_site_rows": rows,
                "assigned": assigned,
                "nearest_d1_per_site": assigned,
            }
        )
    )
    (path.parent / "mapping_and_acceptance.json").write_text(
        json.dumps({"source_sha256": "source"})
    )


def test_condition_comparison_recalculates_positive_gate_and_hashes(
    tmp_path: Path,
) -> None:
    """Available native condition returns actual per-site rows and +1 gate."""
    pilot = tmp_path / "pilot"
    path = (
        pilot
        / "wt_msa_controls_20261127/inspection_export_20261127/boltz_msa_seed20261126/geometry_any_amide.json"
    )
    _write_metric(path)
    result = build_condition_comparison(pilot)
    boltz = next(
        c for c in result["conditions"] if c["case"] == "boltz_msa_seed20261126"
    )
    assert boltz["status"] == "available"
    assert boltz["phase"] == "native"
    assert boltz["positive_phi_gate"]["target_phi_sign"] == 1
    assert boltz["positive_phi_gate"]["input_row_count"] == 48
    assert len(boltz["legacy_assigned_per_site"]) == 4
    assert boltz["source_sha256"]


def test_missing_and_unaccepted_conditions_are_pending(tmp_path: Path) -> None:
    """Missing historical sources remain pending and do not receive fake metrics."""
    result = build_condition_comparison(tmp_path)
    full = next(
        c for c in result["conditions"] if c["condition"] == "full_protein_template"
    )
    masked = next(
        c for c in result["conditions"] if c["condition"] == "corrected_masked_wt"
    )
    assert full["status"] == "pending"
    assert masked["status"] == "pending"
    assert masked["metrics"] is None
    assert result["paired_msa_context"]["status"] == "pending"


def test_rejects_duplicate_or_incomplete_geometry_rows(tmp_path: Path) -> None:
    """Invalid source geometry fails instead of manufacturing a comparison."""
    pilot = tmp_path / "pilot"
    path = (
        pilot
        / "wt_msa_controls_20261127/inspection_export_20261127/boltz_msa_seed20261126/geometry_any_amide.json"
    )
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"all_amide_site_rows": _rows()[:-1], "assigned": []}))
    with pytest.raises(ValueError, match="48 geometry rows"):
        build_condition_comparison(pilot)
