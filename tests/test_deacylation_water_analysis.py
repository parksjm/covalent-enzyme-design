"""Focused geometry and trajectory reader tests for hydration analysis."""

from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pytest
from scipy.io import netcdf_file

from nyl12_covalent.deacylation_water_geometry import (
    _acyl_state,
    _extension_sources,
    _hbond,
    _read_frames,
    _topology,
    angle,
    sha256,
)
from nyl12_covalent.deacylation_water_analysis import analyze


def test_angle_right_angle_and_sha256(tmp_path: Path) -> None:
    """Geometry helper and digest should be deterministic."""
    assert angle(
        np.array([[1.0, 0.0, 0.0]]), np.zeros(3), np.array([[0.0, 1.0, 0.0]])
    ).tolist() == [90.0]
    with pytest.raises(ValueError, match="zero-length"):
        angle(np.zeros((1, 3)), np.zeros(3), np.ones((1, 3)))
    source = tmp_path / "source.bin"
    source.write_bytes(b"abc")
    assert (
        sha256(source)
        == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_water_hbond_geometry_is_minimum_image() -> None:
    """Hydrogen bond helper returns the closest H-A distance and angle."""
    box = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    oxygen = np.array([0.0, 0.0, 0.0])
    hydrogens = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    acceptor = np.array([0.0, 3.0, 0.0])
    ha, da, dha, atom_index = _hbond(oxygen, hydrogens, acceptor[None, :], [42], box)
    assert da == pytest.approx(3.0)
    assert ha == pytest.approx(2.0)
    assert dha == pytest.approx(180.0)
    assert atom_index == 42


def test_hbond_checks_both_equivalent_asp_oxygens() -> None:
    """Closest acceptor geometry can select OD2 over OD1."""
    box = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    oxygen = np.array([0.0, 0.0, 0.0])
    hydrogens = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    acceptors = np.array([[4.0, 2.0, 0.0], [0.0, 3.0, 0.0]])
    ha, da, _, atom_index = _hbond(oxygen, hydrogens, acceptors, [14559, 14560], box)
    assert da == pytest.approx(3.0)
    assert ha == pytest.approx(2.0)
    assert atom_index == 14560


def test_read_frames_preserves_time_stage_and_box(tmp_path: Path) -> None:
    """NetCDF frame reader keeps source time and periodic-cell values."""
    path = tmp_path / "trajectory.nc"
    with netcdf_file(path, "w") as nc:
        nc.createDimension("frame", 1)
        nc.createDimension("atom", 1)
        nc.createDimension("spatial", 3)
        nc.createDimension("cell", 3)
        nc.createVariable("coordinates", "f", ("frame", "atom", "spatial"))[:] = [
            [[1, 2, 3]]
        ]
        nc.createVariable("time", "f", ("frame",))[:] = [2.5]
        nc.createVariable("cell_lengths", "d", ("frame", "cell"))[:] = [[20, 21, 22]]
        nc.createVariable("cell_angles", "d", ("frame", "cell"))[:] = [[90, 90, 90]]
    frame = _read_frames(path, "equil")[0]
    assert frame["stage"] == "equil"
    assert frame["time_ps"] == 2.5
    assert frame["xyz"].tolist() == [[1.0, 2.0, 3.0]]
    assert frame["box"].tolist() == [20.0, 21.0, 22.0, 90.0, 90.0, 90.0]


def test_topology_water_mapping() -> None:
    """Topology parser identifies three-site waters with one-based residues."""
    path = Path(
        "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance/inputs/system.prmtop"
    )
    topology, waters = _topology(path)
    assert len(topology.atoms) == 104922
    assert len(waters) == 28245
    assert waters[0]["residue"] == 1537
    assert waters[0]["o"] == 20187
    assert len(waters[0]["h"]) == 2


def test_accepted_extension_requires_every_scientific_gate(tmp_path: Path) -> None:
    """Accepted extensions with incomplete gate records fail explicitly."""
    segment = tmp_path / "extend_01"
    segment.mkdir()
    (segment / "status.json").write_text('{"status":"accepted"}')
    (segment / "scientific_acceptance.json").write_text('{"numeric": {}}')
    with pytest.raises(KeyError):
        _extension_sources(tmp_path)


def test_extension_collection_keeps_all_nine_and_rejects_accepted_gaps(
    tmp_path: Path,
) -> None:
    """Accepted sampling is uncapped; missing accepted predecessors are errors."""
    record = dict(
        numeric=dict(
            completed_steps=10000,
            final_temperature_K=300,
            electronic_history=dict(verified_scc_history=True),
        ),
        boundary=dict(frozen_max_displacement_A=0, box_unchanged=True),
        C_OG_A=1.34,
        acyl_retained=True,
        finite_temperature_passed=True,
        native_process_passed=True,
        scc_verified=True,
    )
    for index in range(1, 10):
        segment = tmp_path / f"extend_{index:02d}"
        segment.mkdir()
        (segment / "status.json").write_text(json.dumps(dict(status="accepted")))
        (segment / "scientific_acceptance.json").write_text(json.dumps(record))
        (segment / "trajectory.nc").write_bytes(b"collector checks presence only")
    future = tmp_path / "extend_100"
    future.mkdir()
    (future / "status.json").write_text(json.dumps(dict(status="prepared")))
    assert [stage for _, stage in _extension_sources(tmp_path)] == [
        f"extend_{index:02d}" for index in range(1, 10)
    ]
    (tmp_path / "extend_04/status.json").write_text(json.dumps(dict(status="running")))
    with pytest.raises(ValueError, match="missing segment"):
        _extension_sources(tmp_path)


def test_acyl_state_accepts_flat_and_historical_nested_schemas() -> None:
    """Both documented producer schemas resolve the same explicit chemistry."""
    flat = {"C_OG_A": 1.34, "acyl_retained": True}
    nested = {"chemistry": {"C_OG_A": 1.34, "acyl_retained": True}}
    assert _acyl_state(flat) == (1.34, True)
    assert _acyl_state(nested) == (1.34, True)
    assert _acyl_state({**flat, **nested}) == (1.34, True)


@pytest.mark.parametrize(
    "record,error",
    [
        ({"C_OG_A": 1.34}, KeyError),
        ({"C_OG_A": "1.34", "acyl_retained": True}, TypeError),
        ({"C_OG_A": 1.34, "acyl_retained": 1}, TypeError),
        (
            {
                "C_OG_A": 1.34,
                "acyl_retained": True,
                "chemistry": {"C_OG_A": 1.40, "acyl_retained": True},
            },
            ValueError,
        ),
    ],
)
def test_acyl_state_rejects_incomplete_or_conflicting_fields(
    record: dict[str, object], error: type[Exception]
) -> None:
    """Malformed or contradictory scientific chemistry fails explicitly."""
    with pytest.raises(error):
        _acyl_state(record)


def test_full_analysis_matches_endpoint_and_tracks_all_waters(tmp_path: Path) -> None:
    """Production inputs yield finite all-water minima and stage occupancy."""
    base = Path(
        "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
    )
    output = tmp_path / "analysis"
    summary = analyze(output, base)
    expected_stage_counts = {"endpoint": 1, "heat": 20, "equil_300K": 100}
    accepted_paths = []
    for directory in sorted(base.glob("extend_*")):
        status = directory / "status.json"
        if (
            status.is_file()
            and json.loads(status.read_text()).get("status") == "accepted"
        ):
            trajectory = directory / "trajectory.nc"
            with netcdf_file(trajectory, "r", mmap=False) as stream:
                expected_stage_counts[directory.name] = stream.variables[
                    "coordinates"
                ].shape[0]
            accepted_paths.append(trajectory)
    assert summary["frame_count_by_stage"] == expected_stage_counts
    assert summary["frame_count_including_endpoint"] == sum(
        expected_stage_counts.values()
    )
    assert summary["top5_waters"][0]["min_distance_time_ps"] >= 0.0
    for water in summary["top5_waters"]:
        assert (
            0 <= water["near_attack_frame_count"] <= sum(expected_stage_counts.values())
        )
        assert set(water["near_attack_fraction_by_stage"]) == set(expected_stage_counts)
        assert all(
            0 <= fraction <= 1
            for fraction in water["near_attack_fraction_by_stage"].values()
        )
    import csv

    with (output / "nearest_water_per_frame.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert float(rows[0]["O_C_A"]) == pytest.approx(6.245937, abs=1e-5)
    assert len(rows) == sum(expected_stage_counts.values())
    distances = np.array([float(row["O_C_A"]) for row in rows])
    assert np.isfinite(distances).all() and np.all(distances > 0)
    assert min(distances) == pytest.approx(
        summary["top5_waters"][0]["min_attack_distance_A"]
    )
    with netcdf_file(accepted_paths[-1], "r", mmap=False) as stream:
        assert float(rows[-1]["time_ps"]) == pytest.approx(
            float(stream.variables["time"][-1])
        )
    with (output / "per_frame_metrics.csv").open() as stream:
        details = csv.DictReader(stream)
        assert {
            "C_OG1_A",
            "C_Nleaving_A",
            "waterO_Asp272_D_A",
            "waterH_Asp272_A_A",
            "waterO_H_Asp272_angle_deg",
            "waterH_Asp272_closest_acceptor_atom_index_1based",
        }.issubset(details.fieldnames or [])
