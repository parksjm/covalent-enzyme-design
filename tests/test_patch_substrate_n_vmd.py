"""Verify presentation patches preserve bundle integrity and are idempotent."""

import hashlib
from pathlib import Path

import pytest

from scripts.export_direct_pt_progress import write_checksum_manifest
from scripts.patch_substrate_n_vmd import LABEL, patch_bundle


def test_patch_preserves_other_labels_checksums_and_originals(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    vmd = bundle / "reaction_path.vmd"
    original = "label add Bonds $molid/13904 $molid/13916\nlabel textsize 1.0\n"
    vmd.write_text(original)
    (bundle / "reaction_path.dcd").write_bytes(b"unchanged trajectory")
    write_checksum_manifest(bundle)
    old_manifest = (bundle / "MANIFEST.sha256").read_bytes()
    backup = tmp_path / "original"
    journal = tmp_path / "actions.jsonl"
    assert patch_bundle(bundle, backup, journal)
    assert vmd.read_text() == original.replace(
        "label textsize", LABEL + "\nlabel textsize"
    )
    assert (backup / vmd.name).read_text() == original
    assert (backup / "MANIFEST.sha256").read_bytes() == old_manifest
    for line in (bundle / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ")
        assert hashlib.sha256((bundle / name).read_bytes()).hexdigest() == digest
    assert not patch_bundle(bundle, backup, journal)
    assert len(journal.read_text().splitlines()) == 2


@pytest.mark.parametrize("corrupt", ["manifest", "anchor"])
def test_patch_refuses_inconsistent_bundle(tmp_path: Path, corrupt: str) -> None:
    vmd = tmp_path / "reaction_path.vmd"
    vmd.write_text("label textsize 1.0\n" if corrupt == "manifest" else "unexpected\n")
    write_checksum_manifest(tmp_path)
    if corrupt == "manifest":
        vmd.write_text(vmd.read_text() + "changed\n")
    with pytest.raises(ValueError):
        patch_bundle(tmp_path, tmp_path / "backup", tmp_path / "actions.jsonl")
    assert not (tmp_path / "backup").exists()
