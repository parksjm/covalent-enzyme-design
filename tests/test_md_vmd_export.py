"""Tests for compact neutral-replica VMD exports."""

from pathlib import Path

import pytest

from nyl12_covalent.md_vmd_export import (
    cpptraj_topology_input,
    cpptraj_trajectory_input,
    export_replica,
    vmd_script,
)


def test_cpptraj_input_orders_segments_images_and_strips(tmp_path: Path) -> None:
    text = cpptraj_trajectory_input(
        tmp_path / "system.prmtop",
        [tmp_path / "10ns.nc", tmp_path / "90ns.nc"],
        [1000, 9000],
        25,
        tmp_path / "movie.dcd",
    )
    assert text.index("10ns.nc") < text.index("90ns.nc")
    assert "trajin " + str((tmp_path / "90ns.nc").resolve()) + " 9000 9000" in text
    assert "autoimage anchor :1313-1316" in text
    assert "strip :WAT,Na+,Cl-,NA,CL" in text


def test_cpptraj_input_rejects_invalid_inventory(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="inventories"):
        cpptraj_trajectory_input(
            tmp_path / "top", [tmp_path / "one.nc"], [1, 2], 1, tmp_path / "out"
        )
    with pytest.raises(ValueError, match="positive"):
        cpptraj_trajectory_input(
            tmp_path / "top", [tmp_path / "one.nc"], [1], 0, tmp_path / "out"
        )


def test_topology_and_vmd_scripts_match_display_contract(tmp_path: Path) -> None:
    topology = cpptraj_topology_input(tmp_path / "full.parm7", tmp_path / "stripped.parm7")
    assert "parmstrip :WAT,Na+,Cl-,NA,CL" in topology
    viewer = vmd_script(401)
    assert "Expected 401 frames" in viewer
    assert 'mol selection "resid 1313 to 1316 and element C"' in viewer
    assert "mol color ColorID 3" in viewer
    assert "mol showrep $molid 0 0" in viewer
    assert "mol selupdate 3 $molid on" in viewer


def test_vmd_script_requires_frames() -> None:
    with pytest.raises(ValueError, match="positive"):
        vmd_script(0)


def test_export_replica_writes_verified_portable_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise manifest, viewer, provenance, and checksum publication."""
    replica = tmp_path / "neutral_replica_01"
    replica.mkdir()
    for name in ("system.prmtop", "05_production_10ns.nc", "06_production_90ns.nc"):
        (replica / name).write_bytes(name.encode())
    cpptraj = tmp_path / "cpptraj"
    cpptraj.write_text("executable")
    output = replica / "vmd_stripped_100ns"

    def fake_frame_count(path: Path) -> int:
        return 1000 if path.name.startswith("05_") else 9000

    def fake_run(_executable: Path, input_path: Path, log_path: Path) -> None:
        log_path.write_text("verified\n")
        if input_path.name == "trajectory.cpptraj":
            (output / "production_100ns_stripped.dcd").write_bytes(b"dcd")
        else:
            (output / "stripped.parm7").write_text("topology")

    monkeypatch.setattr(
        "nyl12_covalent.md_vmd_export.netcdf_frame_count", fake_frame_count
    )
    monkeypatch.setattr("nyl12_covalent.md_vmd_export._run_cpptraj", fake_run)
    monkeypatch.setattr(
        "nyl12_covalent.md_vmd_export.dcd_counts", lambda _path: (401, 19967)
    )
    manifest = export_replica(replica, cpptraj, output)
    assert manifest["frames"] == 401
    assert manifest["stripped_atoms"] == 19967
    assert manifest["display_frame_interval_ps"] == 250.0
    assert manifest["exact_final_frame_appended"] is True
    assert (output / "view_production.vmd").is_file()
    assert (output / "MANIFEST.sha256").read_text().count("\n") == 9
    assert "visualization only" in (output / "manifest.json").read_text()


def test_export_replica_rejects_unexpected_frame_inventory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    replica = tmp_path / "neutral_replica_01"
    replica.mkdir()
    for name in ("system.prmtop", "05_production_10ns.nc", "06_production_90ns.nc"):
        (replica / name).write_text(name)
    cpptraj = tmp_path / "cpptraj"
    cpptraj.write_text("executable")
    monkeypatch.setattr(
        "nyl12_covalent.md_vmd_export.netcdf_frame_count", lambda _path: 2
    )
    with pytest.raises(ValueError, match=r"1000 \+ 9000"):
        export_replica(replica, cpptraj, replica / "bundle")
