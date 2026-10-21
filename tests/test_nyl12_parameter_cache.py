"""Tests for reproducible Nyl12 parameter artifact caching."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from nyl12_covalent.parameter_cache import ParameterArtifactCache


def test_identity_order_and_round_trip(tmp_path: Path) -> None:
    """Equivalent mapping order yields one verified, reusable entry."""
    source = tmp_path / "parameters.mol2"
    source.write_text("parameter data\n")
    cache = ParameterArtifactCache(tmp_path / "cache")
    identity = {"graph": "abc", "tool": {"version": "1"}}
    reordered = {"tool": {"version": "1"}, "graph": "abc"}

    key = cache.store(identity, {"ligand/parameters.mol2": source})
    assert key == cache.key_for(reordered)
    target = tmp_path / "materialized"
    paths = cache.materialize(reordered, target)

    assert paths == [target / "ligand/parameters.mol2"]
    assert paths[0].read_text() == "parameter data\n"
    manifest = json.loads(
        (tmp_path / "cache" / "parameter_artifacts" / "manifest.json").read_text()
    )
    assert manifest["identity"] == identity


def test_miss_and_existing_entry_reuse(tmp_path: Path) -> None:
    cache = ParameterArtifactCache(tmp_path / "cache")
    identity = {"state": "neutral"}
    assert cache.materialize(identity, tmp_path / "out") is None
    source = tmp_path / "input"
    source.write_bytes(b"original")
    key = cache.store(identity, {"params.dat": source})
    source.write_bytes(b"changed")
    with pytest.raises(RuntimeError, match="different artifacts"):
        cache.store(identity, {"params.dat": source})
    source.write_bytes(b"original")
    assert cache.store(identity, {"params.dat": source}) == key
    result = cache.materialize(identity, tmp_path / "out")
    assert result is not None and result[0].read_bytes() == b"original"


def test_corruption_fails_loudly(tmp_path: Path) -> None:
    source = tmp_path / "input"
    source.write_text("original")
    cache = ParameterArtifactCache(tmp_path / "cache")
    identity = {"charge_method": "example"}
    cache.store(identity, {"params.dat": source})
    (tmp_path / "cache" / "parameter_artifacts" / "params.dat").write_text("tampered")
    with pytest.raises(RuntimeError, match="corrupt"):
        cache.materialize(identity, tmp_path / "out")


def test_collisions_and_unsafe_names_are_rejected(tmp_path: Path) -> None:
    source = tmp_path / "input"
    source.write_text("original")
    cache = ParameterArtifactCache(tmp_path / "cache")
    identity = {"conformer": 3}
    cache.store(identity, {"params.dat": source})
    target = tmp_path / "out"
    target.mkdir()
    (target / "params.dat").write_text("keep")
    with pytest.raises(FileExistsError):
        cache.materialize(identity, target)
    with pytest.raises(ValueError, match="Unsafe"):
        cache.store({"bad": True}, {"../escape": source})


def test_non_json_identity_and_empty_artifacts_fail(tmp_path: Path) -> None:
    with pytest.raises(TypeError):
        ParameterArtifactCache.key_for({"unsupported": object()})
    cache = ParameterArtifactCache(tmp_path / "cache")
    with pytest.raises(ValueError, match="At least one"):
        cache.store({"valid": True}, {})


def test_json_sequence_identity_is_canonicalized(tmp_path: Path) -> None:
    source = tmp_path / "input"
    source.write_text("artifact")
    cache = ParameterArtifactCache(tmp_path / "cache")
    key = cache.store({"coordinates": (1, 2)}, {"params.dat": source})
    assert cache.materialize({"coordinates": [1, 2]}, tmp_path / "out") is not None
    assert key == cache.key_for({"coordinates": [1, 2]})


def test_labeled_entry_preserves_hash_identity_and_reuses_cache(
    tmp_path: Path,
) -> None:
    source = tmp_path / "input"
    source.write_text("artifact")
    cache = ParameterArtifactCache(tmp_path / "cache")
    identity = {"state": "neutral", "conformer": 1}

    key = cache.store(identity, {"params.dat": source}, label="neutral_thr_conf_01")

    entry = tmp_path / "cache" / "neutral_thr_conf_01"
    manifest = json.loads((entry / "manifest.json").read_text())
    assert key == cache.key_for(identity)
    assert manifest["key"] == key
    assert cache.materialize(identity, tmp_path / "out") is not None
    assert cache.store(identity, {"params.dat": source}, label="another_label") == key
    assert not (tmp_path / "cache" / "another_label").exists()


def test_hash_named_legacy_entry_is_reused_and_can_be_migrated(
    tmp_path: Path,
) -> None:
    source = tmp_path / "input"
    source.write_text("legacy artifact")
    cache = ParameterArtifactCache(tmp_path / "cache")
    identity = {"state": "legacy"}
    key = cache.store(identity, {"params.dat": source}, label="temporary_label")
    labeled = tmp_path / "cache" / "temporary_label"
    legacy = tmp_path / "cache" / key
    os.replace(labeled, legacy)
    manifest_before = (legacy / "manifest.json").read_bytes()

    assert cache.store(identity, {"params.dat": source}, label="new_label") == key
    assert legacy.is_dir()
    migrated = cache.migrate_entry(identity, "readable_legacy_entry")

    assert migrated == tmp_path / "cache" / "readable_legacy_entry"
    assert not legacy.exists()
    assert (migrated / "manifest.json").read_bytes() == manifest_before
    record = json.loads(
        (tmp_path / "cache" / ".migrations" / f"{key}.json").read_text()
    )
    assert record["old_directory"] == key
    assert record["new_directory"] == "readable_legacy_entry"
    assert record["key"] == key
    assert cache.materialize(identity, tmp_path / "after_migration") is not None


def test_labels_are_safe_and_collisions_fail(tmp_path: Path) -> None:
    cache = ParameterArtifactCache(tmp_path / "cache")
    first = tmp_path / "first"
    first.write_text("first")
    second = tmp_path / "second"
    second.write_text("second")

    cache.store({"identity": 1}, {"params.dat": first}, label="stable_name")
    with pytest.raises(FileExistsError, match="already in use"):
        cache.store({"identity": 2}, {"params.dat": second}, label="stable_name")
    with pytest.raises(ValueError, match="Unsafe"):
        cache.store({"identity": 3}, {"params.dat": second}, label="../escape")


def test_duplicate_identity_entries_fail_closed(tmp_path: Path) -> None:
    source = tmp_path / "input"
    source.write_text("artifact")
    cache = ParameterArtifactCache(tmp_path / "cache")
    identity = {"state": "duplicate"}
    key = cache.store(identity, {"params.dat": source}, label="first_label")
    duplicate = tmp_path / "cache" / "second_label"
    os.replace(tmp_path / "cache" / "first_label", duplicate)
    (tmp_path / "cache" / key).mkdir()
    (tmp_path / "cache" / key / "manifest.json").write_bytes(
        (duplicate / "manifest.json").read_bytes()
    )
    (tmp_path / "cache" / key / "params.dat").write_bytes(source.read_bytes())

    with pytest.raises(RuntimeError, match="Multiple parameter cache entries"):
        cache.materialize(identity, tmp_path / "out")
