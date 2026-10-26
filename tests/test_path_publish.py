"""Tests for safely refreshing the small downloaded inspection package."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from nyl12_covalent.path_publish import publish_download

ARTIFACTS = {
    "reaction_path.dcd": "reaction_path.dcd",
    "reaction_path.vmd": "reaction_path.vmd",
    "stripped.parm7.gz": "stripped.parm7.gz",
    "qm_indices.tcl": "qm_indices.tcl",
    "nyl12_es_vmd_common.tcl": "nyl12_es_vmd_common.tcl",
    "nyl12_reactive_distances.tcl": "nyl12_reactive_distances.tcl",
    "reaction_distances.png": "distances.png",
    "reaction_distances.json": "distances.json",
    "reaction_distances.csv": "distances.csv",
    "reaction_energies.png": "energies.png",
    "reaction_energies.json": "energies.json",
    "reaction_energies.csv": "energies.csv",
}


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _package(tmp_path: Path) -> tuple[Path, Path, bytes]:
    generation = tmp_path / "generation"
    destination = tmp_path / "download"
    generation.mkdir()
    destination.mkdir()
    for source_name, destination_name in ARTIFACTS.items():
        (generation / source_name).write_bytes(f"new {source_name}".encode())
        (destination / destination_name).write_bytes(f"old {destination_name}".encode())
    original_manifest = (
        "".join(
            f"{_sha256((destination / name).read_bytes())}  {name}\n"
            for name in ARTIFACTS.values()
        )
        + f"{_sha256(b'keep me')}  unrelated.txt\n"
    ).encode()
    (destination / "unrelated.txt").write_bytes(b"keep me")
    (destination / "MANIFEST.sha256").write_bytes(original_manifest)
    return generation, destination, original_manifest


def test_publish_download_preserves_other_entries_and_replaces_atomically(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    generation, destination, old_manifest = _package(tmp_path)
    real_replace = Path.replace
    replacements: list[str] = []

    def observe_replace(source: Path, target: Path) -> Path:
        replacements.append(source.name)
        if source.name == ".MANIFEST.pending":
            for source_name, destination_name in ARTIFACTS.items():
                assert (destination / destination_name).read_bytes() == (
                    generation / source_name
                ).read_bytes()
        else:
            assert (
                destination / target.name
            ).read_bytes() == f"old {target.name}".encode()
        return real_replace(source, target)

    monkeypatch.setattr(Path, "replace", observe_replace)
    publish_download(generation, destination)

    assert (
        generation / "previous_download_manifest.sha256"
    ).read_bytes() == old_manifest
    assert (destination / "unrelated.txt").read_bytes() == b"keep me"
    manifest_lines = (destination / "MANIFEST.sha256").read_text().splitlines()
    manifest = {
        line.split(maxsplit=1)[1]: line.split(maxsplit=1)[0] for line in manifest_lines
    }
    assert manifest["unrelated.txt"] == _sha256(b"keep me")
    for source_name, destination_name in ARTIFACTS.items():
        expected = (generation / source_name).read_bytes()
        assert (destination / destination_name).read_bytes() == expected
        assert manifest[destination_name] == _sha256(expected)
    assert replacements == [f".{name}.pending" for name in ARTIFACTS.values()] + [
        ".MANIFEST.pending"
    ]
    assert not list(destination.glob(".*.pending"))


def test_publish_download_missing_source_leaves_checksum_manifest_unchanged(
    tmp_path: Path,
) -> None:
    generation, destination, old_manifest = _package(tmp_path)
    (generation / "reaction_distances.json").unlink()

    with pytest.raises(FileNotFoundError):
        publish_download(generation, destination)

    assert (destination / "MANIFEST.sha256").read_bytes() == old_manifest
    assert (
        generation / "previous_download_manifest.sha256"
    ).read_bytes() == old_manifest


def test_publish_download_rejects_corrupt_staged_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    generation, destination, old_manifest = _package(tmp_path)

    def corrupt_copy(source: Path, target: Path) -> None:
        target.write_bytes(b"corrupt staged bytes")

    monkeypatch.setattr("nyl12_covalent.path_publish.shutil.copy2", corrupt_copy)
    with pytest.raises(ValueError, match="digest mismatch"):
        publish_download(generation, destination)

    assert (destination / "MANIFEST.sha256").read_bytes() == old_manifest
    assert (destination / "reaction_path.dcd").read_bytes() == b"old reaction_path.dcd"


def test_publish_download_initializes_new_destination(tmp_path: Path) -> None:
    generation = tmp_path / "generation"
    destination = tmp_path / "new_download"
    generation.mkdir()
    for name in ARTIFACTS:
        (generation / name).write_bytes(f"new {name}".encode())

    publish_download(generation, destination)

    assert (generation / "previous_download_manifest.sha256").read_text() == ""
    manifest = (destination / "MANIFEST.sha256").read_text()
    for source_name, destination_name in ARTIFACTS.items():
        assert (destination / destination_name).read_bytes() == (
            generation / source_name
        ).read_bytes()
        assert destination_name in manifest
