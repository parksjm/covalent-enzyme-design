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
    assert (destination / "vmd/view_path.vmd").read_bytes() == (
        latest / "reaction_path.vmd"
    ).read_bytes()
    assert (destination / "plots/distances.png").read_bytes() == (
        latest / "reaction_distances.png"
    ).read_bytes()
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


def test_live_child_has_its_own_checksum_boundary(tmp_path: Path) -> None:
    """A separately published relaxation must not invalidate the parent manifest."""
    config_path, destination, _latest = _fixture(tmp_path)
    config = json.loads(config_path.read_text())
    config["independent_packages"] = ["phase1_unrestrained"]
    config_path.write_text(json.dumps(config))
    child = destination / "phase1_unrestrained"
    child.mkdir(parents=True)
    movie = child / "optimization.dcd"
    movie.write_bytes(b"first live snapshot")
    build_qmmm_package(config_path)
    manifest = (destination / "MANIFEST.sha256").read_text()
    assert "phase1_unrestrained/" not in manifest
    movie.write_bytes(b"new live snapshot")
    for line in manifest.splitlines():
        expected, relative = line.split("  ", 1)
        assert _hash((destination / relative).read_bytes()) == expected


def test_independent_package_cannot_escape_root(tmp_path: Path) -> None:
    """Reject ambiguous or escaping checksum boundaries."""
    config_path, _destination, _latest = _fixture(tmp_path)
    config = json.loads(config_path.read_text())
    config["independent_packages"] = ["../other"]
    config_path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="child directories"):
        build_qmmm_package(config_path)


def test_build_qmmm_package_rejects_unmapped_stage(tmp_path: Path) -> None:
    config_path, _destination, _latest = _fixture(tmp_path, stage="stage2")

    with pytest.raises(KeyError, match="stage2"):
        build_qmmm_package(config_path)


def test_changed_label_configuration_refreshes_existing_package(tmp_path: Path) -> None:
    config_path, destination, latest = _fixture(tmp_path)
    build_qmmm_package(config_path)
    config = json.loads(config_path.read_text())
    config["measurement_atoms"] = {
        "O": {"index": 1, "name": "OG1"},
        "H": {"index": 2, "name": "HG1"},
    }
    config["stage_measurements"] = {"stage1": [["O", "H"]]}
    helper = tmp_path / "labels.tcl"
    helper.write_text("# native reactive labels")
    config["source_files"]["vmd/nyl12_reactive_distances.tcl"] = str(helper)
    config_path.write_text(json.dumps(config))
    build_qmmm_package(config_path)
    viewer = destination / "vmd/view_path.vmd"
    assert "::nyl12_reactive::install" in viewer.read_text()
    assert "{1 OG1 2 HG1}" in viewer.read_text()
    assert (
        destination / "vmd/nyl12_reactive_distances.tcl"
    ).read_bytes() == helper.read_bytes()
    source = json.loads((destination / "package_sources.json").read_text())["sources"][
        "vmd/view_path.vmd"
    ]
    assert source["source_sha256"] == _hash((latest / "reaction_path.vmd").read_bytes())
    assert source["sha256"] == _hash(viewer.read_bytes())


def test_energy_plot_is_in_plots_and_data_are_downloadable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from nyl12_covalent import path_energies, qmmm_inspection_package

    config_path, destination, _ = _fixture(tmp_path)
    config = json.loads(config_path.read_text())
    config["energy_plot"] = {"python": "explicit-python", "script": "energy_plot.py"}
    config_path.write_text(json.dumps(config))
    rows = [{"window": "stage1_1", "relative_total_kcal_mol": 0.0, "source_sha256": {}}]
    monkeypatch.setattr(
        path_energies, "relay_endpoint_energies", lambda run, endpoints: rows
    )

    def render(command: list[str], **kwargs: object) -> None:
        assert command[:2] == ["explicit-python", "energy_plot.py"]
        assert kwargs["check"] is True
        Path(command[-1]).write_bytes(b"energy image")

    monkeypatch.setattr(qmmm_inspection_package.subprocess, "run", render)
    build_qmmm_package(config_path)
    assert (destination / "plots/energies.png").read_bytes() == b"energy image"
    assert json.loads((destination / "data/energies.json").read_text()) == rows
    assert "relative_total_kcal_mol" in (destination / "data/energies.csv").read_text()
    assert "plots/energies.png" in (destination / "MANIFEST.sha256").read_text()
