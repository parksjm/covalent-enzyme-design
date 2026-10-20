"""Unit tests for the extraction inventory without scientific artifacts."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tarfile
from pathlib import Path

import pytest

from scripts import inventory_covalent_pipeline as inventory_module


def _fake_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Create a small Git tree and substitute its explicit inventory list."""
    root = tmp_path / "repo"
    (root / "scripts").mkdir(parents=True)
    (root / "scripts" / "inventory_covalent_pipeline.py").write_text("generator")
    (root / "src" / "pkg").mkdir(parents=True)
    (root / "src" / "pkg" / "module.py").write_bytes(b"source bytes\n")
    (root / "tests").mkdir()
    (root / "tests" / "test_module.py").write_bytes(b"test bytes\n")
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    subprocess.run(["git", "-C", str(root), "config", "user.email", "test@example.org"], check=True)
    subprocess.run(["git", "-C", str(root), "config", "user.name", "Inventory Test"], check=True)
    (root / "src" / "pkg" / "module.py").write_bytes(b"committed bytes\n")
    subprocess.run(["git", "-C", str(root), "add", "src/pkg/module.py"], check=True)
    subprocess.run(["git", "-C", str(root), "commit", "-qm", "fixture"], check=True)
    (root / "src" / "pkg" / "module.py").write_bytes(b"modified source bytes\n")
    monkeypatch.setattr(inventory_module, "__file__", str(root / "scripts/inventory_covalent_pipeline.py"))
    monkeypatch.setattr(inventory_module, "REQUIRED", {
        "src/pkg/module.py": "source module", "tests/test_module.py": "test module",
    })
    return root


def test_inventory_hashes_bytes_status_and_archive(tmp_path: Path,
                                                    monkeypatch: pytest.MonkeyPatch) -> None:
    """Metadata hashes and archive payload describe identical source bytes."""
    root = _fake_repo(tmp_path, monkeypatch)
    document, files = inventory_module.build_inventory(root)
    indexed = {entry["path"]: entry for entry in document["files"]}
    assert indexed["src/pkg/module.py"]["sha256"] == hashlib.sha256(files["src/pkg/module.py"]).hexdigest()
    assert indexed["src/pkg/module.py"]["git_status"] == "tracked_modified"
    assert indexed["tests/test_module.py"]["git_status"] == "untracked"
    assert document["head_captures_untracked_changes"] is False
    snapshot = tmp_path / "snapshot.tar.gz"
    inventory_module.write_snapshot(snapshot, files)
    with tarfile.open(snapshot, "r:gz") as archive:
        assert {member.name: archive.extractfile(member).read() for member in archive.getmembers()} == files


def test_required_path_symlink_and_overwrite_rejected(tmp_path: Path,
                                                      monkeypatch: pytest.MonkeyPatch) -> None:
    """Missing or escaping paths fail, and CLI refuses existing outputs."""
    root = _fake_repo(tmp_path, monkeypatch)
    (root / "tests/test_module.py").unlink()
    with pytest.raises(FileNotFoundError, match="Required"):
        inventory_module.build_inventory(root)
    (root / "tests/test_module.py").symlink_to(tmp_path / "outside")
    (tmp_path / "outside").write_text("outside")
    with pytest.raises(ValueError, match="regular file"):
        inventory_module.build_inventory(root)
    output = tmp_path / "existing.json"
    output.write_text(json.dumps({"keep": True}))
    with pytest.raises(SystemExit):
        inventory_module.main(["--root", str(root), "--output", str(output)])
    assert json.loads(output.read_text()) == {"keep": True}
