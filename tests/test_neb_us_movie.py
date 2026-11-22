"""Focused validation of native NEB-US movie export planning."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from scipy.io import netcdf_file

from nyl12_covalent import neb_us_movie
from nyl12_covalent.neb_us_movie import (
    _load_json,
    _sample_frames,
    _sha256,
    _validate_block,
    export_part,
)


def test_load_json_rejects_missing_path(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        _load_json(tmp_path / "missing.json")


def test_load_json_rejects_non_object(tmp_path: Path) -> None:
    path = tmp_path / "array.json"
    path.write_text("[]")
    with pytest.raises(ValueError, match="Expected JSON object"):
        _load_json(path)


def test_sample_frames_match_native_five_of_one_hundred() -> None:
    assert _sample_frames(100, 5) == [20, 40, 60, 80, 100]


def test_sample_frames_is_unique_for_uneven_counts() -> None:
    assert _sample_frames(11, 4) == [3, 6, 9, 11]


@pytest.mark.parametrize("frame_count,requested", [(0, 1), (4, 5), (5, 0)])
def test_sample_frames_rejects_invalid_counts(frame_count: int, requested: int) -> None:
    with pytest.raises(ValueError):
        _sample_frames(frame_count, requested)


def test_validate_block_rejects_unfinished_process(tmp_path: Path) -> None:
    block = tmp_path / "block"
    block.mkdir()
    (block / "process.json").write_text(json.dumps({"exit_status": 2}))
    (block / "acceptance.json").write_text("{}")
    with pytest.raises(ValueError, match="process did not succeed"):
        _validate_block(block, "x", 1)


def test_validate_block_requires_valid_acceptance(tmp_path: Path) -> None:
    block = tmp_path / "block"
    block.mkdir()
    (block / "process.json").write_text(json.dumps({"exit_status": 0}))
    (block / "acceptance.json").write_text(json.dumps({"completed_steps": 1}))
    with pytest.raises(ValueError, match="acceptance failed"):
        _validate_block(block, "x", 1)


def test_validate_block_rejects_bad_trajectory_digest(tmp_path: Path) -> None:
    block = tmp_path / "block"
    (block / "observables").mkdir(parents=True)
    (block / "trajectory.nc").write_bytes(b"trajectory")
    (block / "process.json").write_text(json.dumps({"exit_status": 0}))
    accepted = {
        "completed_steps": 10000,
        "final_temperature_K": 300.0,
        "temperature_gate_K": [250, 350],
        "finite_energy_temperature": True,
        "SCC_error_detected": False,
        "frozen_outer_max_displacement_A": 0.0,
    }
    (block / "acceptance.json").write_text(json.dumps(accepted))
    obs = {"source_sha256": {str(block / "trajectory.nc"): "bad"}}
    (block / "observables" / "observables.json").write_text(json.dumps(obs))
    with pytest.raises(ValueError, match="Trajectory hash mismatch"):
        _validate_block(block, "topology", 1)


def test_validate_block_accepts_hashed_native_records(tmp_path: Path) -> None:
    block = tmp_path / "block"
    (block / "observables").mkdir(parents=True)
    trajectory = block / "trajectory.nc"
    with netcdf_file(trajectory, "w") as nc:
        nc.createDimension("frame", 100)
        nc.createDimension("atom", 3)
        nc.createDimension("spatial", 3)
        nc.createVariable("coordinates", "f", ("frame", "atom", "spatial"))
        times = nc.createVariable("time", "f", ("frame",))
        times.units = b"picosecond"
        times[:] = list(range(100))
    topology = tmp_path / "system.prmtop"
    topology.write_text("topology")
    process_path = block / "process.json"
    process_path.write_text(json.dumps({"exit_status": 0}))
    acceptance_path = block / "acceptance.json"
    acceptance_path.write_text(
        json.dumps(
            {
                "completed_steps": 10000,
                "final_temperature_K": 300.0,
                "temperature_gate_K": [250, 350],
                "finite_energy_temperature": True,
                "SCC_error_detected": False,
                "frozen_outer_max_displacement_A": 0.0,
            }
        )
    )
    obs_path = block / "observables" / "observables.json"
    obs_path.write_text(
        json.dumps(
            {
                "frames": 100,
                "source_sha256": {
                    str(trajectory): _sha256(trajectory),
                    str(topology): _sha256(topology),
                },
            }
        )
    )

    source, metadata = _validate_block(block, _sha256(topology), 3)

    assert source == trajectory
    assert len(metadata["times"]) == 100
    assert metadata["trajectory_sha256"] == _sha256(trajectory)
    assert set(metadata["control_sha256"]) == {
        str(process_path),
        str(acceptance_path),
        str(obs_path),
    }


def test_hash_test_fixture_is_stable() -> None:
    assert hashlib.sha256(b"x").hexdigest() == (
        "2d711642b726b04401627ca9fbac32f5c8530fb1903cc4db02258717921a4881"
    )


def test_export_part_orders_windows_and_preserves_qm_water(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    topology = campaign / "system.prmtop"
    topology.write_text("input topology")
    windows = [
        {
            "label": label,
            "bias": {"q_acyl_center": -1.7, "q_h_center": -1.4},
            "stages": [
                {"kind": "production", "label": f"block_{i}"} for i in range(1, 7)
            ],
        }
        for label in ("window_001", "window_002")
    ]
    (campaign / "prepared.json").write_text(
        json.dumps({"topology": str(topology), "windows": windows})
    )
    cpptraj = tmp_path / "cpptraj"
    cpptraj.write_text("stub")
    output = tmp_path / "out"
    monkeypatch.setattr(neb_us_movie, "_sha256", lambda _: "a" * 64)
    monkeypatch.setattr(
        neb_us_movie,
        "_validate_block",
        lambda block, _, __: (
            block / "trajectory.nc",
            {
                "observables": {"qm_atom_ids": [2]},
                "times": list(range(100)),
                "time_units": "ps",
                "trajectory_sha256": "b" * 64,
                "control_sha256": {},
            },
        ),
    )

    def topology_factory(path: str) -> SimpleNamespace:
        atoms = [
            SimpleNamespace(idx=0, name="CA", residue=SimpleNamespace(name="ALA")),
            SimpleNamespace(idx=1, name="O", residue=SimpleNamespace(name="WAT")),
        ]
        if not path.endswith("part.parm7"):
            atoms.append(
                SimpleNamespace(idx=2, name="H1", residue=SimpleNamespace(name="WAT"))
            )
        return SimpleNamespace(atoms=atoms)

    monkeypatch.setattr(neb_us_movie, "AmberParm", topology_factory)

    def fake_cpptraj(script: Path, binary: Path, log: Path) -> None:
        del binary
        log.write_text("success")
        (script.parent / "part.parm7").write_text("topology")
        (script.parent / "part.dcd").write_text("trajectory")

    monkeypatch.setattr(neb_us_movie, "_cpptraj", fake_cpptraj)
    monkeypatch.setattr(neb_us_movie, "dcd_counts", lambda _: (60, 2))
    manifest = export_part(campaign, ["window_002", "window_001"], cpptraj, output)

    assert manifest["windows"] == ["window_001", "window_002"]
    assert manifest["frames"] == 60
    assert manifest["original_atom_ids_1based"] == [1, 2]
    assert manifest["qm_atom_ids_1based"] == [2]
    assert manifest["frame_plan"][0]["source_frame1based"] == 20
    assert manifest["frame_plan"][-1]["output_frame0based"] == 59
    assert (output / "manifest.json").is_file()
    assert "!@2" in (output / "write_trajectory.in").read_text()
