"""Focused checks for sparse-coordinate augmentation of direct-US outputs."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.io import netcdf_file

from nyl12_covalent import direct_us_coordinates as module
from nyl12_covalent.direct_us_coordinates import (
    OLD_DISTANCE_COLUMNS,
    augment_arrays,
    augment_stage,
)
from nyl12_covalent.neb_umbrella_observables import _distances


def _fixture() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict[str, str]]]:
    """Build one triclinic frame and matching original observable row."""
    coordinates = np.array(
        [
            [
                [9.6, 1.0, 1.0],
                [0.8, 1.1, 1.0],
                [9.4, 1.2, 1.0],
                [9.9, 1.2, 1.0],
                [10.2, 1.3, 1.0],
            ]
        ],
        dtype=np.float64,
    )
    times = np.array([np.float32(6.0500002)], dtype=np.float32)
    boxes = np.array([[10.0, 11.0, 12.0, 90.0, 90.0, 60.0]])
    frame = coordinates[0]
    measured = {
        name: float(_distances(frame, [first], [second], boxes[0])[0, 0])
        for name, (first, second) in OLD_DISTANCE_COLUMNS.items()
    }
    rows = [
        {
            "time_ps": str(float(times[0])),
            "q_h": str(measured["r_og1_hg1"] - measured["r_nleave_hg1"]),
            **{key: str(value) for key, value in measured.items()},
        }
    ]
    return coordinates, times, boxes, rows


def test_augment_arrays_triclinic_frame_and_q_partition() -> None:
    coordinates, times, boxes, rows = _fixture()

    augmented = augment_arrays(coordinates, times, boxes, rows)

    assert len(augmented) == 1
    frame = coordinates[0]
    expected = float(_distances(frame, [5], [4], boxes[0])[0, 0])
    assert augmented[0]["r_nalpha_hg1"] == pytest.approx(expected)
    assert augmented[0]["q_act"] + augmented[0]["q_del"] == pytest.approx(
        float(rows[0]["q_h"])
    )


@pytest.mark.parametrize("mismatch", ["time", "distance"])
def test_augment_arrays_rejects_csv_mismatch(mismatch: str) -> None:
    coordinates, times, boxes, rows = _fixture()
    if mismatch == "time":
        rows[0]["time_ps"] = "6.2"
    else:
        rows[0]["r_c_n"] = "9.0"

    with pytest.raises(ValueError, match="does not match CSV"):
        augment_arrays(coordinates, times, boxes, rows)


@pytest.mark.parametrize("hash_kind", ["trajectory", "csv"])
def test_augment_stage_rejects_source_hash_mismatch(
    tmp_path: Path, hash_kind: str
) -> None:
    stage = tmp_path / "stage"
    observable_dir = stage / "observables"
    observable_dir.mkdir(parents=True)
    trajectory = stage / "trajectory.nc"
    trajectory.write_bytes(b"trajectory bytes")
    csv_path = observable_dir / "observables.csv"
    csv_path.write_text("time_ps,q_h,r_c_n,r_c_og1,r_og1_hg1,r_nleave_hg1\n")
    (observable_dir / "observables.json").write_text(
        json.dumps(
            {
                "frames": 100,
                "frames_censored": 0,
                "csv_sha256": (
                    "0" * 64
                    if hash_kind == "csv"
                    else hashlib.sha256(csv_path.read_bytes()).hexdigest()
                ),
                "source_sha256": {
                    str(trajectory.resolve()): (
                        "0" * 64
                        if hash_kind == "trajectory"
                        else hashlib.sha256(trajectory.read_bytes()).hexdigest()
                    )
                },
                "roles_1based": {
                    "carbon": 19916,
                    "leaving_n": 19929,
                    "og1": 13916,
                    "hg1": 13917,
                    "nalpha": 13905,
                },
            }
        )
    )
    (stage / "process.json").write_text('{"exit_status": 0}')
    (stage / "acceptance.json").write_text(
        '{"completed_steps": 10000, "finite_energy_temperature": true, '
        '"SCC_error_detected": false}'
    )

    expected_error = "trajectory hash" if hash_kind == "trajectory" else "CSV hash"
    with pytest.raises(ValueError, match=expected_error):
        augment_stage(stage, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_augment_stage_writes_separate_contract(tmp_path: Path, monkeypatch) -> None:
    stage = tmp_path / "stage"
    observable_dir = stage / "observables"
    observable_dir.mkdir(parents=True)
    monkeypatch.setattr(module, "ATOM_IDS", (1, 2, 3, 4, 5))
    coordinates, _, boxes, _ = _fixture()
    sparse = np.repeat(coordinates, 100, axis=0).astype(np.float32)
    times = np.arange(100, dtype=np.float32) * 0.05 + np.float32(6.05)
    box_lengths = np.repeat(boxes[:, :3], 100, axis=0).astype(np.float32)
    box_angles = np.repeat(boxes[:, 3:], 100, axis=0).astype(np.float32)
    box_values = np.column_stack([box_lengths, box_angles])
    base_row = {
        name: str(float(_distances(sparse[0], [first], [second], box_values[0])[0, 0]))
        for name, (first, second) in OLD_DISTANCE_COLUMNS.items()
    }
    base_row["q_h"] = str(
        float(base_row["r_og1_hg1"]) - float(base_row["r_nleave_hg1"])
    )
    with netcdf_file(str(stage / "trajectory.nc"), "w") as nc:
        nc.createDimension("frame", 100)
        nc.createDimension("atom", 5)
        nc.createDimension("spatial", 3)
        nc.createDimension("cell_spatial", 3)
        nc.createDimension("cell_angular", 3)
        coords_var = nc.createVariable("coordinates", "f", ("frame", "atom", "spatial"))
        time_var = nc.createVariable("time", "f", ("frame",))
        length_var = nc.createVariable("cell_lengths", "f", ("frame", "cell_spatial"))
        angle_var = nc.createVariable("cell_angles", "f", ("frame", "cell_angular"))
        coords_var[:] = sparse
        time_var[:] = times
        length_var[:] = box_lengths
        angle_var[:] = box_angles

    rows = [{"time_ps": str(float(time)), **base_row} for time in times]
    csv_path = observable_dir / "observables.csv"
    columns = ["time_ps", "q_h", *OLD_DISTANCE_COLUMNS]
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    process_path = stage / "process.json"
    acceptance_path = stage / "acceptance.json"
    process_path.write_text('{"exit_status": 0}')
    acceptance_path.write_text(
        '{"completed_steps": 10000, "finite_energy_temperature": true, '
        '"SCC_error_detected": false}'
    )
    trajectory = stage / "trajectory.nc"
    receipt_path = observable_dir / "observables.json"
    receipt_path.write_text(
        json.dumps(
            {
                "frames": 100,
                "frames_censored": 0,
                "csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                "source_sha256": {
                    str(trajectory.resolve()): hashlib.sha256(
                        trajectory.read_bytes()
                    ).hexdigest()
                },
                "roles_1based": {
                    "carbon": 1,
                    "leaving_n": 2,
                    "og1": 3,
                    "hg1": 4,
                    "nalpha": 5,
                },
            }
        )
    )

    output = tmp_path / "augmented"
    receipt = augment_stage(stage, output)

    assert (output / "observables.csv").is_file()
    assert (output / "augmentation.json").is_file()
    assert receipt["rows_count"] == 100
    assert receipt["frames_censored"] == 0
    assert receipt["parity"]["all_four_distances_passed"] is True
    assert receipt["parity"]["q_h_partition_passed"] is True
    assert str(process_path.resolve()) in receipt["source_sha256"]
