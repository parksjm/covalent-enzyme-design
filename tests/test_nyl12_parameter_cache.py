"""Tests for reproducible Nyl12 parameter artifact caching."""

from __future__ import annotations

import json
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
    manifest = json.loads((tmp_path / "cache" / key / "manifest.json").read_text())
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
    key = cache.store(identity, {"params.dat": source})
    (tmp_path / "cache" / key / "params.dat").write_text("tampered")
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
