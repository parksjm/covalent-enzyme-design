"""Check all-site reporting distinguishes failed cases from rejected poses."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from scripts import summarize_original_matched_geometry as target
from nyl12_covalent.negative_amide_gate import evaluate_signed_amide_gate


def rows() -> list[dict[str, Any]]:
    return [
        {
            "site": site,
            "ligand": ligand,
            "carbonyl_c": f"C{i}",
            "carbonyl_o": f"O{i}",
            "leaving_n": f"N{i}",
            "d1": 3.0 + i / 10,
            "d2": 4.0,
            "d3": 3.2,
            "a_bd": 105.0,
            "herschlag_phi": -80.0,
            "tier": "A",
        }
        for site in target.SITES
        for ligand in "IJKL"
        for i in range(3)
    ]


@pytest.fixture
def cohort(tmp_path: Path) -> Path:
    (tmp_path / "controller").mkdir()
    post = tmp_path / "relaxation/good/postprocess"
    post.mkdir(parents=True)
    records = rows()
    # Three sites compete for one eligible ligand; fourth has no eligible amide.
    for row in records:
        if row["ligand"] == "I" and row["site"] != "DH" and row["carbonyl_c"] == "C1":
            row["herschlag_phi"] = 80.0
    for snapshot in ("before", "after"):
        gate = {
            "case": "good",
            "snapshot": snapshot,
            "all_amide_site_rows": records,
            **evaluate_signed_amide_gate(records, target_phi_sign=1),
        }
        (post / f"pose_gate_{snapshot}.json").write_text(json.dumps(gate))
    cases = {
        "good": {"status": "completed", "postprocess": str(post)},
        "bad_native": {"status": "skipped_native_failure"},
        "bad_assembly": {"status": "failed_assembly_precheck"},
    }
    (tmp_path / "controller/summary.json").write_text(json.dumps({"cases": cases}))
    return tmp_path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def test_all_four_sites_full_geometry_and_missing_mm(cohort: Path) -> None:
    summary = target.summarize_cohort(cohort, cohort / "comparison")
    assert summary["measured_mm_cases"] == 1
    assert summary["cases_without_mm_result"] == 2
    assert summary["measured_site_snapshot_rows"] == 8
    assert summary["all_amide_geometry_rows"] == 96
    assert summary["cases"]["good"]["positive_ab_after"] == 1
    assert "positive_ab_after" not in summary["cases"]["bad_assembly"]
    site = read_csv(cohort / "comparison/active_site_comparison.csv")
    assert len(site) == 24
    measured = [r for r in site if r["case"] == "good" and r["snapshot"] == "after"]
    assert {r["site"] for r in measured} == set(target.SITES)
    assert sum(r["positive_ab_assigned"] == "True" for r in measured) == 1
    assert (
        sum(
            r["assignment_status"] == "unassigned_one_to_one_conflict" for r in measured
        )
        == 2
    )
    dh = next(r for r in measured if r["site"] == "DH")
    assert dh["assignment_status"] == "reject_no_passing_amide"
    assert dh["selected_ligand"] == ""
    assert float(dh["diagnostic_nearest_d1_herschlag_phi"]) < 0
    assert all(r["positive_ab_assigned"] == "" for r in site if r["case"] != "good")
    full = read_csv(cohort / "comparison/all_amide_geometry.csv")
    assert len(full) == 96
    assert {r["case"] for r in full} == {"good"}
    assert len(summary["input_sha256"]) == 3
    assert "no fabricated geometry" in (cohort / "comparison/README.md").read_text()


@pytest.mark.parametrize(
    "fault", ["sign", "identity", "snapshot", "assignment", "rows"]
)
def test_bad_gate_fails_before_publication(cohort: Path, fault: str) -> None:
    path = cohort / "relaxation/good/postprocess/pose_gate_before.json"
    gate = json.loads(path.read_text())
    if fault == "sign":
        gate["target_phi_sign"] = -1
    elif fault == "identity":
        gate["case"] = "wrong"
    elif fault == "snapshot":
        gate["snapshot"] = "after"
    elif fault == "assignment":
        gate["assigned_ab_count"] = 3
    else:
        gate["all_amide_site_rows"].pop()
    path.write_text(json.dumps(gate))
    with pytest.raises(ValueError):
        target.summarize_cohort(cohort, cohort / "comparison")
    assert not (cohort / "comparison").exists()


def test_existing_output_preserved(cohort: Path) -> None:
    out = cohort / "comparison"
    out.mkdir()
    with pytest.raises(FileExistsError):
        target.summarize_cohort(cohort, out)


def test_failed_only_has_headers_and_no_fake_geometry(cohort: Path) -> None:
    summary_path = cohort / "controller/summary.json"
    summary_path.write_text(json.dumps({"cases": {"bad": {"status": "failed"}}}))
    summary = target.summarize_cohort(cohort, cohort / "comparison")
    assert summary["measured_site_snapshot_rows"] == 0
    assert read_csv(cohort / "comparison/all_amide_geometry.csv") == []
    assert len(read_csv(cohort / "comparison/active_site_comparison.csv")) == 8


def test_invalid_case_name_fails(cohort: Path) -> None:
    (cohort / "controller/summary.json").write_text(
        json.dumps({"cases": {"../bad": {"status": "failed"}}})
    )
    with pytest.raises(ValueError, match="identity"):
        target.summarize_cohort(cohort, cohort / "comparison")
