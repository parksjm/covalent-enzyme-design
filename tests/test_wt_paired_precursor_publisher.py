"""Tests for the paired-precursor two-WT publication gate."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from design.nyl12_relay_ts_ensemble_20261126.broader_pilot.wt_paired_precursor_controls_20261127 import (
    publish_campaign as publisher,
)


def _campaign() -> tuple[dict, dict, dict, dict, dict]:
    cases = [{"case": case} for case in publisher.EXPECTED_CASES]
    campaign = {
        "cases": cases,
        "prediction_budget": 2,
        "target_phi_sign": 1,
        "template_policy": "none; stock native entry; no structural input",
        "paired_precursor_msa": {"source_sha256": "source-hash"},
    }
    cpu = {
        "status": "accepted",
        "prediction_budget": 2,
        "models": [
            {
                "case": case,
                "record_templates_empty": True,
                "actual_native_template_features_equal_stock_dummy": True,
                "native_nonquery_coherent_paired_rows": 5,
            }
            for case in publisher.EXPECTED_CASES
        ],
    }
    native = {
        "accepted": {
            case: {"status": "native_identity_graph_canonical_roundtrip_accepted"}
            for case in publisher.EXPECTED_CASES
        }
    }
    relax = {
        "completed": {
            case: {"force_converged": True} for case in publisher.EXPECTED_CASES
        },
        "protocol": "synthetic matched protocol",
    }
    target = {"status": "PI_APPROVED_POSITIVE_POSE_REFERENCE", "target_phi_sign": 1}
    return campaign, cpu, native, relax, target


def _rows() -> list[dict]:
    sites = ("AE", "BF", "CG", "DH")
    ligands = "IJKL"
    return [
        {
            "site": sites[i // 12],
            "ligand": ligands[(i // 3) % 4],
            "carbonyl_c": f"C{i}",
            "monomer_boundary": [i % 3 + 1, i % 3 + 2],
            "hypothetical_fragment_lengths": [i % 3 + 1, 3 - i % 3],
            "register_is_product_hypothesis": True,
            "d1": float(i),
            "herschlag_phi": 20.0,
            "tier": "reject",
        }
        for i in range(48)
    ]


def test_select_cases_requires_exactly_two_accepted_controls() -> None:
    """Only the two accepted paired WT predictions and relaxations qualify."""
    campaign, cpu, native, relax, _ = _campaign()
    assert (
        publisher.select_cases(campaign, cpu, native, relax) == publisher.EXPECTED_CASES
    )
    relax["completed"].pop(publisher.EXPECTED_CASES[0])
    with pytest.raises(ValueError, match="Both paired WT Amber relaxations"):
        publisher.select_cases(campaign, cpu, native, relax)


def test_cpu_acceptance_rejects_templates_dummy_mismatch_or_no_pairs() -> None:
    """Native template-free acceptance is explicit for both cases."""
    campaign, cpu, native, relax, _ = _campaign()
    cpu["models"][0]["actual_native_template_features_equal_stock_dummy"] = False
    with pytest.raises(ValueError, match="stock dummy"):
        publisher.select_cases(campaign, cpu, native, relax)
    campaign, cpu, native, relax, _ = _campaign()
    cpu["models"][1]["native_nonquery_coherent_paired_rows"] = 0
    with pytest.raises(ValueError, match="genuine paired homolog"):
        publisher.select_cases(campaign, cpu, native, relax)


def test_pairing_receipt_binds_retained_origins_to_source(tmp_path: Path) -> None:
    """A source hash and unique origin record count are required."""
    source = tmp_path / "source.a3m"
    source.write_text(">query\nAC\n>pair1\nAC\n")
    import hashlib

    receipt = {
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "paired_homolog_row_count": 1,
        "source_row_count": 2,
        "row_origins": [{"key": 1, "source_row": 1, "header": "pair1"}],
    }
    publisher.validate_pairing_receipt(receipt, source)
    receipt["source_sha256"] = "wrong"
    with pytest.raises(ValueError, match="source hash"):
        publisher.validate_pairing_receipt(receipt, source)


def test_live_paired_source_receipt_row_origins_contract() -> None:
    """The staged campaign receipt binds every retained homolog to the source."""
    source = publisher.BASE / "sources/precursor_source.a3m"
    receipt_path = publisher.BASE / "msas/precursor_pairing_receipt.json"
    if not source.is_file() or not receipt_path.is_file():
        pytest.skip("Live precursor MSA source has not been staged")
    receipt = json.loads(receipt_path.read_text())
    publisher.validate_pairing_receipt(receipt, source)
    assert len(receipt["row_origins"]) == receipt["paired_homolog_row_count"]
    assert (
        len(receipt["excluded_rows"])
        == receipt["source_row_count"] - len(receipt["row_origins"]) - 1
    )


def test_zero_and_partial_target_assignment_slots_are_counted_correctly() -> None:
    """Four slots may be empty; assigned count equals populated slots only."""
    case = publisher.EXPECTED_CASES[0]
    gate = {
        "case": case,
        "snapshot": "before",
        "target_phi_sign": 1,
        "input_row_count": 48,
        "all_amide_site_rows": [{"x": i} for i in range(48)],
        "assigned_ab_count": 0,
        "assigned": [None, None, None, None],
    }
    publisher.validate_gate(gate, case, "before")
    gate.update(assigned_ab_count=1, assigned=[{"site": "AE"}, None, None, None])
    publisher.validate_gate(gate, case, "before")
    gate["assigned_ab_count"] = 0
    with pytest.raises(ValueError, match="Assigned-row count"):
        publisher.validate_gate(gate, case, "before")


def test_positive_phi_gate_requires_all48_and_assigned_count() -> None:
    """Reject wrong-sign or truncated before/after geometry gates."""
    case = publisher.EXPECTED_CASES[0]
    gate = {
        "case": case,
        "snapshot": "after",
        "target_phi_sign": 1,
        "input_row_count": 48,
        "all_amide_site_rows": _rows(),
        "assigned_ab_count": 0,
        "assigned": [None, None, None, None],
    }
    publisher.validate_gate(gate, case, "after")
    gate["target_phi_sign"] = -1
    with pytest.raises(ValueError, match="target_phi_sign"):
        publisher.validate_gate(gate, case, "after")


def _write_case(base: Path, case: str) -> None:
    source = base / "native_models" / case
    canonical = source / "canonical"
    canonical.mkdir(parents=True)
    (canonical / "raw_prediction.cif").write_text("data_canonical_pair_input\n#\n")
    (canonical / "mapping_and_acceptance.json").write_text("{}\n")
    (source / "original_native_prediction.cif").write_text(
        "data_original_prediction\n#\n"
    )
    rows = _rows()
    (source / "geometry_all48_and_negative_gate.json").write_text(
        json.dumps(
            {
                "all_amide_site_rows": rows,
                "legacy_assigned": [],
                "nearest_d1_per_site": rows[:4],
                "target_phi_ab_gate": {
                    "target_phi_sign": 1,
                    "input_row_count": 48,
                    "assigned_ab_count": 0,
                    "assigned": [None, None, None, None],
                },
            }
        )
    )
    (base / "inputs").mkdir(exist_ok=True)
    (base / "inputs" / f"{case}.yaml").write_text(
        json.dumps({"version": 1, "sequences": [{"protein": {"sequence": "AC"}}]})
    )
    post = base / "amber_relaxation" / case / "postprocess"
    for phase in ("before", "after"):
        (post / "canonical_before_after" / f"{case}_{phase}").mkdir(parents=True)
        gate = {
            "case": case,
            "snapshot": phase,
            "target_phi_sign": 1,
            "input_row_count": 48,
            "all_amide_site_rows": rows,
            "assigned_ab_count": 0,
            "assigned": [None, None, None, None],
        }
        (post / f"pose_gate_{phase}.json").write_text(json.dumps(gate))
    (post / "case_postprocessing.json").write_text("{}\n")
    (post / "canonical_before_after").mkdir(exist_ok=True)
    (post / "canonical_before_after/before_after_all48.csv").write_text("site\n")


def test_main_publishes_288_rows_and_honest_control_context(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The full publication path writes a portable two-case bundle in tmp_path."""
    base = tmp_path / "paired_campaign"
    base.mkdir()
    monkeypatch.setattr(publisher, "BASE", base)
    monkeypatch.setattr(publisher, "PILOT", tmp_path / "pilot")
    campaign, cpu, native, relax, target = _campaign()
    (base / "campaign.json").write_text(json.dumps(campaign))
    (base / "cpu_acceptance.json").write_text(json.dumps(cpu))
    (base / "collector_state.json").write_text(json.dumps(native))
    (base / "relaxation_state.json").write_text(json.dumps(relax))
    (base / "target_pose_convention.json").write_text(json.dumps(target))
    for case in publisher.EXPECTED_CASES:
        _write_case(base, case)
    (base / "sources").mkdir()
    source = base / "sources/precursor_source.a3m"
    source.write_text(">query\nAC\n>homolog1\nAC\n")
    import hashlib

    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    (base / "sources/source_provenance.json").write_text(
        json.dumps({"status": "accepted", "precursor_a3m_sha256": source_hash})
    )
    (base / "msas").mkdir()
    pairing = {
        "source_sha256": source_hash,
        "source_row_count": 2,
        "paired_homolog_row_count": 1,
        "row_origins": [{"key": 1, "source_row": 1, "header": "homolog1"}],
        "excluded_rows": [],
    }
    (base / "msas/precursor_pairing_receipt.json").write_text(json.dumps(pairing))
    campaign["paired_precursor_msa"]["source_sha256"] = source_hash
    (base / "campaign.json").write_text(json.dumps(campaign))
    pilot = tmp_path / "pilot"
    unpaired = pilot / "wt_msa_controls_20261127/msas"
    unpaired.mkdir(parents=True)
    (unpaired / "alpha_homolog_unpaired.a3m").write_text(">q\nAC\n>h\nAC\n")
    (unpaired / "beta_homolog_unpaired.a3m").write_text(">q\nAC\n")
    template = (
        pilot / "wt_protein_template_controls_20261127/inspection_export_20261127"
    )
    template.mkdir(parents=True)
    (template / "crystal_reference.pdb").write_text(
        "ATOM      1  CA  ALA A  20       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    (template / "MANIFEST.sha256").write_text("")
    (
        pilot / "wt_msa_controls_20261127/calibration_20261127/msa_canonical_roundtrip"
    ).mkdir(parents=True)

    publisher.main()

    output = base / "inspection_export_20261127"
    rows = json.loads((output / "comparison_all_rows.json").read_text())
    assert len(rows) == 288
    assert {row["case"] for row in rows} == set(publisher.EXPECTED_CASES)
    context = json.loads((output / "msa_depth_context.json").read_text())
    assert context["interpretation"].startswith("MSA sources and depths differ")
    references = json.loads((output / "comparison_sources.json").read_text())
    assert references["corrected_masked_wt"]["status"] == "pending"
    assert (output / "MANIFEST.sha256").is_file()
    assert (
        json.loads((base / "publication_acceptance.json").read_text())[
            "comparison_rows"
        ]
        == 288
    )
