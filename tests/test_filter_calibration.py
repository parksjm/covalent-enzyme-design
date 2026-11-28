"""Focused checks for retrospective geometry-filter calibration."""

from copy import deepcopy

import pytest

from nyl12_covalent.filter_calibration import analyze_cohort


def fixture_data(variants=("WT", "F101E"), count=2):
    config = {
        "id": "demo",
        "system": "PA6",
        "register": "r1",
        "experimental_comparison": "substrate-X",
        "variants": list(variants),
        "expected_models": count,
        "thresholds": {
            "d1_max": 3.5,
            "oxyanion_max": 4.0,
            "a1_min": 85,
            "a1_max": 125,
            "herschlag_phi_min": 60,
            "herschlag_phi_max": 125,
        },
    }
    records = []
    for variant in variants:
        for index in range(count):
            records.append(
                {
                    "variant": variant,
                    "model_index": index,
                    "geometry_tier": "B",
                    "d1": 3.5,
                    "d2": 4.0,
                    "d3": 4.0,
                    "a1": 85,
                    "herschlag_phi": 60 if index == 0 else -125,
                    "interface_delta": 1.0 + index,
                }
            )
    return config, records


def test_summary_normalizes_wt_and_keeps_endpoint_labels_separate():
    config, records = fixture_data()
    activities = [
        {
            "system": "PA6",
            "comparison": "substrate-X",
            "variant": "F101E",
            "relative_activity": 0.0,
        },
        {
            "system": "PA6",
            "comparison": "substrate-Y",
            "variant": "F101E",
            "relative_activity": 0.2,
        },
        {
            "system": "PA66",
            "comparison": "substrate-X",
            "variant": "F101E",
            "relative_activity": 9.0,
        },
    ]
    result = analyze_cohort(records, activities, config)
    wt, mutant = result["per_variant"]
    assert wt["selected_endpoint_relative_activity"] == 1.0
    assert wt["endpoint_label"] == "wt_normalized_reference"
    assert mutant["selected_endpoint_relative_activity"] == 0.0
    assert mutant["endpoint_label"] == "zero_reported_not_confirmed_inactive"
    assert mutant["any_product_reported_positive"] is True
    assert result["summary"]["known_positive_variants_with_retained_model"] == 1
    assert result["summary"]["known_positive_variant_count"] == 1
    assert result["summary"]["descriptive_activity_pearson"] is None


def test_missing_endpoint_and_gate_counts_phi_sign_and_inclusive_bounds():
    config, records = fixture_data()
    records[2]["interface_delta"] = -27.85
    records[3]["interface_delta"] = -29.12
    result = analyze_cohort(records, [], config)
    mutant = result["per_variant"][1]
    assert mutant["endpoint_label"] == "missing"
    assert mutant["n_tier_ab"] == 2
    assert mutant["n_phi_positive"] == 1 and mutant["n_phi_negative"] == 1
    assert mutant["n_phi_zero"] == 0
    assert mutant["gate_fail_counts"] == {
        "d1": 0,
        "d2": 0,
        "d3": 0,
        "a1": 0,
        "abs_phi": 0,
    }
    assert mutant["gate_only_failure_counts"]["abs_phi"] == 0
    assert mutant["median_interface_delta_tier_ab"] == pytest.approx(-28.485)
    assert not result["tier_ab_mismatches"]


def test_isolated_gate_failure_and_stored_tier_mismatch_are_reported():
    config, records = fixture_data(count=1)
    records[1]["geometry_tier"] = "B"
    records[1]["d2"] = 4.01
    result = analyze_cohort(records, [], config)
    mutant = result["per_variant"][1]
    assert mutant["gate_fail_counts"]["d2"] == 1
    assert mutant["gate_only_failure_counts"]["d2"] == 1
    assert result["tier_ab_mismatches"] == [
        {"variant": "F101E", "model_index": 0, "stored_tier": "B", "gate_pass": False}
    ]
    assert result["per_model"][1]["gate_flags"]["abs_phi"]


def test_completeness_duplicate_and_nonfinite_inputs_fail_explicitly():
    config, records = fixture_data(count=2)
    with pytest.raises(ValueError, match="Incomplete models"):
        analyze_cohort(records[:-1], [], config)
    with pytest.raises(ValueError, match="Duplicate model"):
        analyze_cohort(records + [records[0]], [], config)
    bad = deepcopy(records)
    bad[0]["d1"] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        analyze_cohort(bad, [], config)


def test_duplicate_selected_activity_and_negative_activity_fail():
    config, records = fixture_data(count=1)
    activity = {
        "system": "PA6",
        "comparison": "substrate-X",
        "variant": "F101E",
        "relative_activity": 0.1,
    }
    with pytest.raises(ValueError, match="Duplicate selected"):
        analyze_cohort(records, [activity, activity], config)
    bad = dict(activity, relative_activity=-0.1)
    with pytest.raises(ValueError, match="nonnegative"):
        analyze_cohort(records, [bad], config)


def test_distance_must_be_nonnegative_and_tier_must_be_valid():
    config, records = fixture_data(count=1)
    bad = deepcopy(records)
    bad[0]["d3"] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        analyze_cohort(bad, [], config)
    bad = deepcopy(records)
    bad[0]["geometry_tier"] = "D"
    with pytest.raises(ValueError, match="geometry_tier"):
        analyze_cohort(bad, [], config)
