"""Focused tests for modular WT QM/MM snapshot extraction helpers."""

from pathlib import Path

import pytest

from scripts.prepare_modular_wt_qmmm_snapshots import (
    cpptraj_text,
    explicit_native_torsion_selections,
    source_record,
)


def test_cpptraj_text_is_single_frame_and_imaged(tmp_path: Path) -> None:
    topology = tmp_path / "system.prmtop"
    trajectory = tmp_path / "production.nc"
    output = tmp_path / "frame.nc"
    text = cpptraj_text(topology, trajectory, 1171, output)
    assert f"trajin {trajectory.resolve()} 1171 1171 1" in text
    assert "autoimage anchor :1313-1316" in text
    assert f"trajout {output.resolve()} netcdf" in text


def test_cpptraj_text_rejects_nonpositive_frame(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="positive"):
        cpptraj_text(tmp_path / "a", tmp_path / "b", 0, tmp_path / "c")


def test_source_record_hashes_file(tmp_path: Path) -> None:
    path = tmp_path / "source"
    path.write_bytes(b"abc")
    record = source_record(path)
    assert record["bytes"] == 3
    assert record["sha256"] == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )
    assert record["path"] == str(path.resolve())


def test_explicit_native_torsion_selection(tmp_path: Path) -> None:
    geometry = tmp_path / "geometry.csv"
    geometry.write_text(
        "frame,time_ns,d1,d2,d3,attack,ligand_rmsd,projected_phi\n"
        "1,0.0,3.0,3.0,3.0,105,2,-90\n"
        "2,0.05,3.4,3.0,3.0,105,2,-90\n",
        encoding="utf-8",
    )
    torsions = tmp_path / "torsions.csv"
    torsions.write_text(
        "frame,tier,native_torsions,psi_delta,chi1_delta\n"
        "1,A,True,-10,5\n2,B,False,-50,5\n",
        encoding="utf-8",
    )
    result = explicit_native_torsion_selections(geometry, torsions, [1])
    assert result[0]["frame"] == 1 and result[0]["tier"] == "A"
    with pytest.raises(ValueError, match="fails joint"):
        explicit_native_torsion_selections(geometry, torsions, [2])
