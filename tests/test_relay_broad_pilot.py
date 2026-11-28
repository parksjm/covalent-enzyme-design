"""Tests for the bounded broad-pilot preparation helpers."""

from __future__ import annotations

import pytest

from nyl12_covalent.relay_broad_pilot import sample_schedule, symmetry_manifest


def test_symmetry_manifest_ties_all_four_alpha_and_beta_copies() -> None:
    manifest = symmetry_manifest({"alpha": [2, 7], "beta": [3]})

    assert manifest["symmetry_residues"] == (
        "A2,B2,C2,D2|A7,B7,C7,D7|E3,F3,G3,H3"
    )
    assert manifest["symmetry_weights"] == "1,0,0,0|1,0,0,0|1,0,0,0"
    assert set(manifest["redesigned_residues"].split()) == {
        "A2", "B2", "C2", "D2", "A7", "B7", "C7", "D7",
        "E3", "F3", "G3", "H3",
    }


@pytest.mark.parametrize("remodeled_count,expected", [(0, 12), (1, 15), (4, 24)])
def test_schedule_respects_fixed_raw_sample_ceiling(
    remodeled_count: int, expected: int
) -> None:
    originals = ["medoid01", "medoid03"]
    remodeled = [f"rfd3_{index}" for index in range(remodeled_count)]

    schedule = sample_schedule(originals, remodeled)

    assert sum(item["raw_samples"] for item in schedule) == expected
    assert [item["raw_samples"] for item in schedule[:2]] == [6, 6]
    assert all(item["raw_samples"] == 3 for item in schedule[2:])


def test_schedule_rejects_more_than_four_remodeled_parents() -> None:
    with pytest.raises(ValueError, match="Zero to four"):
        sample_schedule(
            ["medoid01", "medoid03"], [f"rfd3_{index}" for index in range(5)]
        )
