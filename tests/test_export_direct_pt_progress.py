"""Tests for direct proton-transfer VMD progress packaging."""

import json
from pathlib import Path

import pytest

from scripts import export_direct_pt_progress as progress_export
from scripts.export_direct_pt_progress import (
    completed_restarts,
    export,
    vmd_script,
    write_checksum_manifest,
)


def test_completed_restarts_follow_recorded_order(tmp_path: Path) -> None:
    run = tmp_path / "run"
    for label in ("pt_01", "pt_02"):
        folder = run / "01_coordinate_drive" / label
        folder.mkdir(parents=True)
        (folder / "optimized.ncrst").write_text(label)
    (run / "status.json").write_text(
        json.dumps(
            {
                "state": "running_coordinate_drive",
                "completed_windows": ["pt_01", "pt_02"],
            }
        )
    )
    assert [path.parent.name for path in completed_restarts(run)] == [
        "pt_01",
        "pt_02",
    ]


def test_attack_only_export_prepends_start_and_exports_completed_windows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = tmp_path / "run"
    start = tmp_path / "starting.ncrst"
    start.write_text("start")
    windows = run / "windows"
    for label in ("attack_01", "attack_02"):
        folder = windows / label
        folder.mkdir(parents=True)
        (folder / "optimized.ncrst").write_text(label)
    release = run / "release" / "optimized.ncrst"
    release.parent.mkdir()
    release.write_text("release")
    (run / "status.json").write_text(
        json.dumps(
            {
                "state": "completed",
                "starting_restart": str(start),
                "completed_windows": ["attack_01", "attack_02"],
            }
        )
    )
    paths = completed_restarts(run, "attack_only")
    assert paths == [
        start,
        windows / "attack_01/optimized.ncrst",
        windows / "attack_02/optimized.ncrst",
    ]

    input_files: list[str] = []

    def fake_run_cpptraj(cpptraj: Path, input_path: Path, cwd: Path, log: Path) -> None:
        input_files.append(input_path.read_text())
        log.write_text("complete")
        if input_path.name == "trajectory.cpptraj":
            (cwd / "reaction_path.dcd").write_text("fake trajectory")
        else:
            (cwd / "stripped.parm7").write_text("fake topology")

    monkeypatch.setattr(progress_export, "run_cpptraj", fake_run_cpptraj)
    monkeypatch.setattr(progress_export, "dcd_counts", lambda _path: (3, 42))
    topology = tmp_path / "topology.parm7"
    topology.write_text("topology")
    cpptraj = tmp_path / "cpptraj"
    cpptraj.write_text("executable placeholder")

    generation = export(run, topology, cpptraj, tmp_path / "exports", "attack_only")

    trajectory_input = input_files[0]
    assert trajectory_input.count("trajin ") == 3
    assert f"trajin {start.resolve()}" in trajectory_input
    assert (
        f"trajin {(windows / 'attack_01/optimized.ncrst').resolve()}"
        in trajectory_input
    )
    assert (
        f"trajin {(windows / 'attack_02/optimized.ncrst').resolve()}"
        in trajectory_input
    )
    assert str(release.resolve()) not in trajectory_input
    manifest = json.loads((generation / "manifest.json").read_text())
    assert manifest["frames"] == 3
    assert manifest["coordinate"] == "attack_only"


def test_asynchronous_restarts_include_unrestrained_endpoint(tmp_path: Path) -> None:
    """The asynchronous movie contains driven endpoints plus the release test."""
    run = tmp_path / "run"
    driven = run / "01_coupled_coordinate_drive/asynchronous_01"
    release = run / "02_unrestrained_ti1_endpoint"
    driven.mkdir(parents=True)
    release.mkdir(parents=True)
    (driven / "optimized.ncrst").write_text("driven")
    (release / "optimized.ncrst").write_text("release")
    (run / "status.json").write_text(
        json.dumps(
            {
                "state": "completed",
                "completed_windows": ["asynchronous_01"],
            }
        )
    )
    paths = completed_restarts(run, "asynchronous")
    assert [path.parent.name for path in paths] == [
        "asynchronous_01",
        "02_unrestrained_ti1_endpoint",
    ]


def test_asynchronous_restarts_prefer_active_extension_over_old_release(
    tmp_path: Path,
) -> None:
    """An extension movie remains continuous and omits the old failed release."""
    run = tmp_path / "run"
    driven = run / "01_coupled_coordinate_drive/asynchronous_24"
    old_release = run / "02_unrestrained_ti1_endpoint"
    extension = run / "03_coordinate_drive_extension_to_ti1/asynchronous_25"
    for directory in (driven, old_release, extension):
        directory.mkdir(parents=True)
        (directory / "optimized.ncrst").write_text(directory.name)
    (run / "status.json").write_text(
        json.dumps({"state": "completed", "completed_windows": ["asynchronous_24"]})
    )
    (run / "extension_status.json").write_text(
        json.dumps(
            {
                "state": "running",
                "completed_windows": ["asynchronous_25"],
            }
        )
    )
    paths = completed_restarts(run, "asynchronous")
    assert [path.parent.name for path in paths] == [
        "asynchronous_24",
        "asynchronous_25",
    ]


def test_completed_restarts_rejects_unknown_coordinate(tmp_path: Path) -> None:
    """Unknown path layouts fail explicitly."""
    with pytest.raises(ValueError, match="Unsupported reaction-path coordinate"):
        completed_restarts(tmp_path, "unknown")


def test_vmd_script_has_required_hidden_and_reactive_representations() -> None:
    text = vmd_script(10)
    assert "mol showrep $molid 0 0" in text
    assert "DynamicBonds" in text
    assert "remove_reactive_bond $molid 13915 13916" in text
    assert "13915 13916 1.40" in text
    assert "13904 13916 1.40" in text
    assert "13915 19915 1.90" in text
    assert "{13915 13916 1.40}" not in text
    assert "resid 1313 to 1316" in text
    assert "mol color Element" in text
    assert "resid 1313 to 1316 and element C and not index" in text
    assert "label add Bonds $molid/13915 $molid/13916" in text
    assert "label add Bonds $molid/19916 $molid/1637" in text
    assert "label add Bonds $molid/19916 $molid/2731" in text
    assert "label add Angles $molid/13915 $molid/19915 $molid/19916" in text
    assert "label add Angles $molid/13916 $molid/13915" not in text
    assert "Loaded 10 finalized pathway frames" in text


def test_vmd_script_rejects_empty_path() -> None:
    with pytest.raises(ValueError, match="positive"):
        vmd_script(0)


def test_vmd_script_accepts_named_umbrella_trajectory() -> None:
    text = vmd_script(80, "umbrella_windows.dcd", "umbrella-sampling")
    assert "umbrella_windows.dcd" in text
    assert "Loaded 80 umbrella-sampling frames" in text


def test_vmd_script_rejects_nonlocal_trajectory_name() -> None:
    with pytest.raises(ValueError, match="local DCD"):
        vmd_script(10, "../frames.dcd")


def test_checksum_manifest_covers_bundle_files(tmp_path: Path) -> None:
    (tmp_path / "a.vmd").write_text("a")
    (tmp_path / "b.dcd").write_text("b")
    write_checksum_manifest(tmp_path)
    text = (tmp_path / "MANIFEST.sha256").read_text()
    assert "  a.vmd\n" in text
    assert "  b.dcd\n" in text
    assert "MANIFEST.sha256" not in text
