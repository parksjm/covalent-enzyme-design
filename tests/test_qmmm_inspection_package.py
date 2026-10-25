"""Tests for the human-readable QM/MM inspection package layout."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from nyl12_covalent.qmmm_inspection_package import build_qmmm_package

PACKAGE_FILES = (
    "reaction_path.dcd",
    "reaction_path.vmd",
    "stripped.parm7.gz",
    "qm_indices.tcl",
    "nyl12_es_vmd_common.tcl",
    "reaction_distances.png",
    "reaction_distances.csv",
    "reaction_distances.json",
)


def _hash(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _fixture(tmp_path: Path, stage: str = "stage1") -> tuple[Path, Path, Path]:
    run = tmp_path / "scientific_run"
    latest = run / "live_inspection" / "latest"
    latest.mkdir(parents=True)
    for name in PACKAGE_FILES:
        (latest / name).write_bytes(f"live {name}".encode())
    data = {
        "identity": "inspection-identity-1",
        "rows": [{"window": "stage1_1", "stage": stage, "accepted": True}],
    }
    (latest / "reaction_distances.json").write_text(json.dumps(data))
    window = run / "stage1_1"
    window.mkdir()
    for name in ("acceptance.json", "process.json", "minimize.mdin", "restraints.RST"):
        (window / name).write_text(f"window {name}")
    source_files = {}
    for target, source_name in (
        ("scripts/run_path.py", "run_path.py"),
        ("scripts/initial_driver.py", "initial_driver.py"),
    ):
        source = tmp_path / source_name
        source.write_text(f"# {source_name}\n")
        source_files[target] = str(source)
    destination = tmp_path / "qmmm_package"
    config = {
        "run": str(run),
        "destination": str(destination),
        "stage_names": {"stage1": "01_og1_to_asp272"},
        "source_files": source_files,
    }
    config_path = tmp_path / "layout.json"
    config_path.write_text(json.dumps(config))
    return config_path, destination, latest


def test_build_qmmm_package_maps_files_hashes_and_is_idempotent(
    tmp_path: Path,
) -> None:
    config_path, destination, latest = _fixture(tmp_path)
    assert build_qmmm_package(config_path) == destination

    window = destination / "01_og1_to_asp272/window_01"
    assert (window / "geometry.json").read_text() == "window acceptance.json"
    assert (window / "run.json").read_text() == "window process.json"
    assert (window / "optimize.mdin").read_text() == "window minimize.mdin"
    assert (window / "restraints.RST").read_text() == "window restraints.RST"
    assert (destination / "vmd/view_path.vmd").read_bytes() == (latest / "reaction_path.vmd").read_bytes()
    assert (destination / "plots/distances.png").read_bytes() == (latest / "reaction_distances.png").read_bytes()
    assert (destination / "scripts/run_path.py").is_file()
    order = json.loads((destination / "window_order.json").read_text())
    assert order == [
        {
            "vmd_frame": 0,
            "directory": "01_og1_to_asp272/window_01",
            "original_window": "stage1_1",
            "accepted_current_policy": True,
        }
    ]
    readme = (destination / "README.md").read_text()
    assert "plots/distances.png" in readme
    assert "01_og1_to_asp272/window_00" in readme
    assert "02_asp272_to_nterm/window_00" in readme

    sources = json.loads((destination / "package_sources.json").read_text())["sources"]
    for relative, metadata in sources.items():
        packaged = destination / relative
        assert metadata["sha256"] == _hash(packaged.read_bytes())
    manifest = {
        line.split(maxsplit=1)[1]: line.split(maxsplit=1)[0]
        for line in (destination / "MANIFEST.sha256").read_text().splitlines()
    }
    for packaged in destination.rglob("*"):
        if packaged.is_file() and packaged.name != "MANIFEST.sha256":
            relative = str(packaged.relative_to(destination))
            assert manifest[relative] == _hash(packaged.read_bytes())

    saved_manifest = (destination / "MANIFEST.sha256").read_bytes()
    saved_manifest_mtime = (destination / "MANIFEST.sha256").stat().st_mtime_ns
    saved_movie = (destination / "vmd/reaction_path.dcd").read_bytes()
    (latest / "reaction_path.dcd").write_bytes(b"changed source after publish")
    assert build_qmmm_package(config_path) == destination
    assert (destination / "MANIFEST.sha256").read_bytes() == saved_manifest
    assert (destination / "MANIFEST.sha256").stat().st_mtime_ns == saved_manifest_mtime
    assert (destination / "vmd/reaction_path.dcd").read_bytes() == saved_movie


def test_build_qmmm_package_fails_for_missing_live_inspection(
    tmp_path: Path,
) -> None:
    config_path, _destination, _latest = _fixture(tmp_path)
    config = json.loads(config_path.read_text())
    config["run"] = str(tmp_path / "missing_run")
    config_path.write_text(json.dumps(config))

    with pytest.raises(FileNotFoundError):
        build_qmmm_package(config_path)


def test_build_qmmm_package_rejects_unmapped_stage(tmp_path: Path) -> None:
    config_path, _destination, _latest = _fixture(tmp_path, stage="stage2")

    with pytest.raises(KeyError, match="stage2"):
        build_qmmm_package(config_path)
