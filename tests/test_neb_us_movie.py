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
        _validate_block(block, "x", 1, [1])


def test_validate_block_requires_valid_acceptance(tmp_path: Path) -> None:
    block = tmp_path / "block"
    block.mkdir()
    (block / "process.json").write_text(json.dumps({"exit_status": 0}))
    (block / "acceptance.json").write_text(json.dumps({"completed_steps": 1}))
    with pytest.raises(ValueError, match="Unknown production acceptance schema"):
        _validate_block(block, "x", 1, [1])


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
        _validate_block(block, "topology", 1, [1])


def test_validate_block_accepts_hashed_native_records(tmp_path: Path) -> None:
    block = tmp_path / "block"
    (block / "observables").mkdir(parents=True)
    trajectory = block / "trajectory.nc"
    with netcdf_file(trajectory, "w") as nc:
        nc.createDimension("frame", 100)
        nc.createDimension("atom", 3)
        nc.createDimension("spatial", 3)
        coordinates = nc.createVariable(
            "coordinates", "f", ("frame", "atom", "spatial")
        )
        coordinates[:] = [[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]] * 100
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

    source, metadata = _validate_block(
        block, _sha256(topology), 3, [20, 40, 60, 80, 100]
    )

    assert source == trajectory
    assert len(metadata["times"]) == 100
    assert metadata["trajectory_sha256"] == _sha256(trajectory)
    assert set(metadata["control_sha256"]) == {
        str(process_path),
        str(acceptance_path),
        str(obs_path),
    }
    acceptance_path.write_text(
        json.dumps(
            {
                "accepted": True,
                "completed_steps": 10000,
                "final_temperature_K": 300.0,
                "temperature_gate_K": [250, 350],
                "finite_energy_temperature": True,
                "every_integration_step_SCC_converged": True,
                "electronic_history": {
                    "passed": True,
                    "verified_scc_history": True,
                    "recovered_trials": [{"trial": 1, "recovered": True}],
                    "unresolved_trials": [],
                    "fatal_markers": [],
                },
                "frozen_outer_max_displacement_A": 0.0,
            }
        )
    )
    _validate_block(block, _sha256(topology), 3, [20])

    with netcdf_file(trajectory, "a") as nc:
        nc.variables["coordinates"][:] = 0.0
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
    with pytest.raises(ValueError, match="coordinates are invalid"):
        _validate_block(block, _sha256(topology), 3, [20])


def test_hash_test_fixture_is_stable() -> None:
    assert hashlib.sha256(b"x").hexdigest() == (
        "2d711642b726b04401627ca9fbac32f5c8530fb1903cc4db02258717921a4881"
    )


@pytest.mark.parametrize("use_continuation", [False, True])
def test_export_part_orders_windows_and_preserves_qm_water(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, use_continuation: bool
) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    topology = campaign / "system.prmtop"
    topology.write_text("input topology")
    windows = [
        {
            "label": label,
            "chart": "activation" if label == "window_001" else "delivery",
            "bias": {
                "feature_names": ["q_acyl", "q_act"],
                "centers_A": [-1.7, -1.4],
            },
            "stages": [
                {"kind": "production", "label": f"block_{i}"} for i in range(1, 7)
            ],
        }
        for label in ("window_001", "window_002")
    ]
    (campaign / "prepared.json").write_text(
        json.dumps({"topology": str(topology), "windows": windows})
    )
    override = None
    if use_continuation:
        override = tmp_path / "continuation_prepared.json"
        override.write_text((campaign / "prepared.json").read_text())
        (campaign / "prepared.json").write_text("{}")
    cpptraj = tmp_path / "cpptraj"
    cpptraj.write_text("stub")
    output = tmp_path / "out"
    monkeypatch.setattr(neb_us_movie, "_sha256", lambda _: "a" * 64)
    monkeypatch.setattr(
        neb_us_movie,
        "_validate_block",
        lambda block, _, __, ___: (
            block / "trajectory.nc",
            {
                "observables": {"qm_atom_ids": [2]},
                "roles_1based": {
                    "carbon": 2,
                    "leaving_n": 2,
                    "og1": 2,
                    "hg1": 2,
                    "nalpha": 2,
                },
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
    monkeypatch.setattr(neb_us_movie, "dcd_counts", lambda _: (40, 2))
    manifest = export_part(
        campaign,
        ["window_002", "window_001"],
        cpptraj,
        output,
        production_blocks=4,
        prepared_manifest=override,
    )

    assert manifest["windows"] == ["window_001", "window_002"]
    assert manifest["window_charts"] == {
        "window_001": "activation",
        "window_002": "delivery",
    }
    assert manifest["frame_plan"][0]["bias_centers_A"] == {
        "q_acyl": -1.7,
        "q_act": -1.4,
    }
    if use_continuation:
        assert str(override) in manifest["source_sha256"]
        assert (campaign / "prepared.json").read_text() == "{}"
    assert manifest["frames"] == 40
    assert manifest["original_atom_ids_1based"] == [1, 2]
    assert manifest["qm_atom_ids_1based"] == [2]
    assert manifest["frame_plan"][0]["source_frame1based"] == 20
    assert manifest["frame_plan"][-1]["output_frame0based"] == 39
    assert (output / "manifest.json").is_file()
    assert "!@2" in (output / "write_trajectory.in").read_text()
