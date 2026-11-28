"""Verify namespace changes preserve every raw lineage through deduplication."""

from copy import deepcopy

import pytest

from nyl12_covalent.relay_ensemble_sequences import pool_samples
from nyl12_covalent.relay_pilot_namespace import namespace_pool


def test_prefix_precedes_deduplication_and_retains_sources() -> None:
    """Duplicate sequences retain both origin IDs and one coherent candidate."""
    raw = [
        dict(sample_id="parent1_1", sequence="AAA"),
        dict(sample_id="parent2_1", sequence="AAA"),
        dict(sample_id="parent2_2", sequence="BBB"),
    ]
    result = namespace_pool(raw, pool=pool_samples, namespace="prot101")
    assert [x["candidate_id"] for x in result] == [
        "relay_prot101_01",
        "relay_prot101_02",
    ]
    assert result[0]["lineage"] == ["prot101_parent1_1", "prot101_parent2_1"]
    assert result[0]["sample_id"] == raw[0]["sample_id"]
    assert result[0]["source_sample_id"] == "parent1_1"
    assert result[0]["source_candidate_id"] == "relay_01"
    assert all(x["sample_id"].startswith("prot101_") for x in raw)
    with pytest.raises(ValueError, match="unnamespaced"):
        namespace_pool(raw, pool=pool_samples, namespace="prot101")


@pytest.mark.parametrize("name", ["Bad", "bad_name", "a space", ""])
def test_invalid_namespace_preserves_input(name: str) -> None:
    """Invalid naming cannot partially modify raw native evidence."""
    raw = [dict(sample_id="parent_1", sequence="AAA")]
    original = deepcopy(raw)
    with pytest.raises(ValueError, match="namespace"):
        namespace_pool(raw, pool=pool_samples, namespace=name)
    assert raw == original


def test_duplicate_raw_ids_fail() -> None:
    """Every independent sample must have a unique pre-pooling identity."""
    with pytest.raises(ValueError, match="unique"):
        namespace_pool(
            [dict(sample_id="same"), dict(sample_id="same")],
            pool=pool_samples,
            namespace="prot101",
        )


def test_missing_or_invented_lineage_fails() -> None:
    """No absent representative or lost raw sample can pass namespace export."""
    raw = [dict(sample_id="a", sequence="AAA")]
    with pytest.raises(ValueError, match="representative"):
        namespace_pool(
            raw,
            pool=lambda s: [
                dict(candidate_id="relay_01", sample_id=s[0]["sample_id"], lineage=[])
            ],
            namespace="prot101",
        )
    raw = [dict(sample_id="a", sequence="AAA")]
    with pytest.raises(ValueError, match="lost"):
        namespace_pool(raw, pool=lambda _: [], namespace="prot101")
