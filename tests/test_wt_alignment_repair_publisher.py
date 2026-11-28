"""Focused tests for feature-qualified six-case publication."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from design.nyl12_relay_ts_ensemble_20261126.broader_pilot.boltz_masked_loop95_designs_20261127.template_alignment_repair01 import (
    publish_corrected_comparison as publisher,
)


def _states() -> tuple[dict, dict]:
    native_ok = "native_identity_graph_canonical_roundtrip_accepted"
    native = {
        "accepted": {case: {"status": native_ok} for case in publisher.REPAIRED_CASES}
    }
    relax = {
        "completed": {
            case: {"force_converged": True} for case in publisher.REPAIRED_CASES
        }
    }
    return native, relax


def _feature_receipt() -> dict:
    return {
        "status": "accepted",
        "expected_query_mask_enforced": True,
        "all_retained_query_author_indices_exact": True,
        "matched_repaired_features_across_six": True,
        "repair_prediction_budget": 6,
    }


def _row(index: int = 0) -> dict:
    return {
        "site": ("AE", "BF", "CG", "DH")[index // 12],
        "ligand": "IJKL"[(index // 3) % 4],
        "carbonyl_c": f"C{index:02d}",
        "monomer_boundary": [index % 3 + 1, index % 3 + 2],
        "hypothetical_fragment_lengths": [index % 3 + 1, 3 - index % 3],
        "register_is_product_hypothesis": True,
    }


def _all_rows() -> list[dict]:
    return [_row(index) for index in range(48)]


def test_select_final_cases_is_exactly_repaired_six() -> None:
    """Selection accepts only the six repaired case IDs."""
    native, relax = _states()
    assert publisher.select_final_cases(native, relax) == publisher.REPAIRED_CASES
    native["accepted"][publisher.REPAIRED_CASES[0]]["status"] = "failed"
    with pytest.raises(ValueError, match="Native acceptance"):
        publisher.select_final_cases(native, relax)


def test_missing_repaired_case_fails_closed() -> None:
    """A missing repaired seed cannot produce a partial comparison."""
    native, relax = _states()
    native["accepted"].pop(publisher.REPAIRED_CASES[-1])
    with pytest.raises(ValueError, match="Required corrected-six inputs missing"):
        publisher.select_final_cases(native, relax)


def test_feature_acceptance_requires_all_six_matched_features() -> None:
    """Wrong feature scope or budget blocks publication."""
    receipt = _feature_receipt()
    publisher.validate_feature_acceptance(receipt)
    receipt["matched_repaired_features_across_six"] = False
    with pytest.raises(ValueError, match="matched_repaired_features_across_six"):
        publisher.validate_feature_acceptance(receipt)
    receipt = _feature_receipt()
    receipt["repair_prediction_budget"] = 2
    with pytest.raises(ValueError, match="repair_prediction_budget"):
        publisher.validate_feature_acceptance(receipt)


def test_positive_phi_gate_requires_complete_register() -> None:
    """Reject wrong sign, wrong case, or truncated register evidence."""
    case = publisher.REPAIRED_CASES[0]
    gate = {
        "case": case,
        "snapshot": "after",
        "target_phi_sign": 1,
        "input_row_count": 48,
        "all_amide_site_rows": [_row()] * 48,
        "assigned_ab_count": 0,
        "assigned": [None] * 4,
    }
    publisher.validate_gate(gate, case, "after")
    gate["target_phi_sign"] = -1
    with pytest.raises(ValueError, match="target_phi_sign"):
        publisher.validate_gate(gate, case, "after")


@pytest.mark.parametrize("assigned_count", range(5))
def test_assigned_ab_count_matches_nonnull_slots(assigned_count: int) -> None:
    """The declared assignment count equals occupied positions among four slots."""
    case = publisher.REPAIRED_CASES[0]
    assigned = [
        ({"index": index} if index < assigned_count else None) for index in range(4)
    ]
    gate = {
        "case": case,
        "snapshot": "after",
        "target_phi_sign": 1,
        "input_row_count": 48,
        "all_amide_site_rows": _all_rows(),
        "assigned_ab_count": assigned_count,
        "assigned": assigned,
    }
    publisher.validate_gate(gate, case, "after")
    gate["assigned_ab_count"] = (assigned_count + 1) % 5
    with pytest.raises(ValueError, match="Assigned-row count mismatch"):
        publisher.validate_gate(gate, case, "after")


def test_assigned_gate_requires_four_slots() -> None:
    """A truncated or expanded assigned slot list is rejected."""
    case = publisher.REPAIRED_CASES[0]
    gate = {
        "case": case,
        "snapshot": "after",
        "target_phi_sign": 1,
        "input_row_count": 48,
        "all_amide_site_rows": _all_rows(),
        "assigned_ab_count": 0,
        "assigned": [None] * 3,
    }
    with pytest.raises(ValueError, match="four assigned slots"):
        publisher.validate_gate(gate, case, "after")


def _write_case(root: Path, case: str) -> None:
    inputs = root / "inputs"
    inputs.mkdir(exist_ok=True)
    input_case = case.removesuffix("_repair01")
    (inputs / f"{input_case}.yaml").write_text(f"case: {case}\n")
    source = root / "native_models" / case
    canonical = source / "canonical"
    canonical.mkdir(parents=True)
    (canonical / "raw_prediction.cif").write_text("data_pair\n#\n")
    rows = _all_rows()
    (canonical / "mapping_and_acceptance.json").write_text(json.dumps({"all48": rows}))
    (source / "original_native_prediction.cif").write_text("data_original\n#\n")
    (source / "geometry_all48_and_negative_gate.json").write_text(
        json.dumps(
            {
                "all_amide_site_rows": rows,
                "target_phi_ab_gate": {"target_phi_sign": 1},
            }
        )
    )
    post = root / "amber_relaxation" / case / "postprocess"
    for phase in ("before", "after"):
        (post / "canonical_before_after" / f"{case}_{phase}").mkdir(parents=True)
        gate = {
            "case": case,
            "snapshot": phase,
            "target_phi_sign": 1,
            "input_row_count": 48,
            "all_amide_site_rows": rows,
            "assigned_ab_count": 0,
            "assigned": [None] * 4,
        }
        (post / f"pose_gate_{phase}.json").write_text(json.dumps(gate))
    (post / "case_postprocessing.json").write_text("{}\n")
    (post / "canonical_before_after").mkdir(exist_ok=True)
    (post / "canonical_before_after/before_after_all48.csv").write_text("case\n")


def test_main_publishes_only_repaired_six_from_temp_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The complete packaging path creates a manifest without original models."""
    repair = tmp_path / "template_alignment_repair01"
    repair.mkdir()
    temp_publisher = repair / "publisher.py"
    temp_publisher.write_text(Path(publisher.__file__).read_text())
    monkeypatch.setattr(publisher, "__file__", str(temp_publisher))
    (repair / "native_feature_acceptance.json").write_text(
        json.dumps(_feature_receipt())
    )
    inventory = {
        "status": "accepted",
        "source_atoms": 8940,
        "invented_coordinates": 0,
        "golden_author_lists_verified_against_actual_source_CIF": True,
    }
    (repair / "source_author_inventory_acceptance.json").write_text(
        json.dumps(inventory)
    )
    coordinate_identity = {
        "status": "accepted",
        "all_source_heavy_atom_coordinates_match_actual_CIF": 8940,
        "all_source_CA_author_coordinates_and_AA_match": 1196,
        "repaired_processed_template_hashes_identical_across_six": True,
    }
    (repair / "source_coordinate_identity_acceptance.json").write_text(
        json.dumps(coordinate_identity)
    )
    native, relax = _states()
    relax["protocol"] = {"method": "synthetic matched protocol fixture"}
    (repair / "collector_state.json").write_text(json.dumps(native))
    (repair / "relaxation_state.json").write_text(json.dumps(relax))
    for case in publisher.REPAIRED_CASES:
        _write_case(repair, case)
    reference = (
        repair.parent.parent
        / "dgx_wt_control_qualification_20261127/inspection_export_20261127"
    )
    reference.mkdir(parents=True)
    (reference / "crystal_reference.pdb").write_text(
        "ATOM      1  CA  ALA A  20       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    (repair.parent / "template_alignment_qualification.json").write_text(
        json.dumps({"limitation": "synthetic exploratory limitation"})
    )
    protection = repair.parent.parent
    (protection / "canonical_protection_hardened_20261127.json").write_text(
        json.dumps({"policy": "fixture"})
    )
    (protection / "protection_requalification_20261127.json").write_text(
        json.dumps({"status": "fixture"})
    )

    publisher.main()

    output = repair / "inspection_export_20261127"
    assert (output / "MANIFEST.sha256").is_file()
    assert (output / "source_author_inventory_acceptance.json").is_file()
    assert (output / "source_coordinate_identity_acceptance.json").is_file()
    assert (output / "canonical_protection_hardened_20261127.json").is_file()
    assert (output / "protection_requalification_20261127.json").is_file()
    assert "diagnostic-only" in (output / "README.md").read_text()
    assert (
        "does not establish that either residue is autocatalytic"
        in (output / "README.md").read_text()
    )
    summary = json.loads((output / "campaign_summary.json").read_text())
    assert summary["cases"] == list(publisher.REPAIRED_CASES)
    assert (
        len(json.loads((output / "comparison_all_rows.json").read_text())) == 6 * 3 * 48
    )
    assert not any("WT_seed_20261126_native" == path.name for path in output.iterdir())
    assert (
        json.loads((repair / "publication_acceptance.json").read_text())[
            "old_wt_in_comparison"
        ]
        is False
    )
