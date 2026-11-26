"""Add sparse native-coordinate relay diagnostics to accepted direct-US stages."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.io import netcdf_file

from nyl12_covalent.neb_umbrella_observables import _distances

ATOM_IDS = (19916, 19929, 13916, 13917, 13905)
ATOM_NAMES = ("carbon", "leaving_n", "og1", "hg1", "nalpha")
OLD_DISTANCE_COLUMNS = {
    "r_c_n": (1, 2),
    "r_c_og1": (1, 3),
    "r_og1_hg1": (3, 4),
    "r_nleave_hg1": (2, 4),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def augment_arrays(
    coordinates: np.ndarray,
    times: np.ndarray,
    boxes: np.ndarray,
    rows: Sequence[Mapping[str, str]],
) -> list[dict[str, float]]:
    """Compute and cross-check direct-US reaction-coordinate diagnostics.

    Coordinates use the sparse atom order carbon, leaving N, OG1, HG1, Nalpha.
    Boxes use Amber lengths/angles. Existing CSV rows are checked before any
    new diagnostic is returned.

    Args:
        coordinates: Sparse coordinates with shape (frames, 5, 3), in Angstrom.
        times: Native time values, shape (frames,).
        boxes: Amber cell lengths and angles, shape (frames, 6).
        rows: Existing CSV rows in frame order.

    Returns:
        Per-frame sparse distances and q_act/q_del values.

    Raises:
        ValueError: If shapes, times, or existing reactive distances disagree.
    """
    xyz = np.asarray(coordinates)
    time_values = np.asarray(times)
    cell_values = np.asarray(boxes)
    if xyz.ndim != 3 or xyz.shape[1:] != (5, 3) or not np.isfinite(xyz).all():
        raise ValueError("Sparse coordinates must be finite with shape (frames, 5, 3)")
    frame_count = xyz.shape[0]
    if (
        frame_count == 0
        or time_values.shape != (frame_count,)
        or cell_values.shape != (frame_count, 6)
        or len(rows) != frame_count
        or not np.isfinite(time_values).all()
        or not np.isfinite(cell_values).all()
    ):
        raise ValueError(
            "Sparse coordinate, time, box, and CSV frame counts must match"
        )

    output: list[dict[str, float]] = []
    for index, (frame, time_value, box, row) in enumerate(
        zip(xyz, time_values, cell_values, rows, strict=True)
    ):
        if abs(float(row["time_ps"]) - float(time_value)) > 1e-6:
            raise ValueError(f"Native time does not match CSV at frame {index + 1}")
        measured: dict[str, float] = {}
        for name, (first, second) in OLD_DISTANCE_COLUMNS.items():
            measured[name] = float(_distances(frame, [first], [second], box)[0, 0])
            if abs(measured[name] - float(row[name])) > 1e-8:
                raise ValueError(
                    f"Native {name} does not match CSV at frame {index + 1}"
                )
        measured["r_nalpha_hg1"] = float(_distances(frame, [5], [4], box)[0, 0])
        measured["q_act"] = measured["r_og1_hg1"] - measured["r_nalpha_hg1"]
        measured["q_del"] = measured["r_nalpha_hg1"] - measured["r_nleave_hg1"]
        if abs(measured["q_act"] + measured["q_del"] - float(row["q_h"])) > 1e-8:
            raise ValueError(
                f"Relay coordinates do not sum to q_h at frame {index + 1}"
            )
        output.append(measured)
    return output


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def augment_stage(stage: Path, output: Path) -> dict[str, Any]:
    """Write coordinate-derived diagnostics for one accepted production stage.

    The native trajectory and existing observables are read-only. The output
    directory must not exist; all input identities are pinned in its receipt.

    Args:
        stage: One ``03_production_block_*_5ps`` directory.
        output: New, separate directory for augmented CSV and provenance.

    Returns:
        Receipt written to ``output/augmentation.json``.

    Raises:
        FileNotFoundError: If a required input is absent.
        ValueError: If acceptance, hashes, native shape, or cross-checks fail.
        FileExistsError: If output already exists.
    """
    stage = Path(stage).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    if output == stage or stage in output.parents:
        raise ValueError("Augmented output must be outside the source stage")
    csv_path = stage / "observables" / "observables.csv"
    receipt_path = stage / "observables" / "observables.json"
    trajectory = stage / "trajectory.nc"
    process_path = stage / "process.json"
    acceptance_path = stage / "acceptance.json"
    receipt = _read_json(receipt_path)
    process = _read_json(process_path)
    acceptance = _read_json(acceptance_path)
    source_hashes = receipt.get("source_sha256", {})
    raw_expected = source_hashes.get(str(trajectory))
    raw_hash = _sha256(trajectory)
    if not raw_expected or raw_hash != raw_expected:
        raise ValueError("Native trajectory hash does not match observables receipt")
    csv_hash = _sha256(csv_path)
    if csv_hash != receipt.get("csv_sha256"):
        raise ValueError("Observable CSV hash does not match observables receipt")
    if (
        process.get("exit_status") != 0
        or acceptance.get("completed_steps") != 10000
        or acceptance.get("finite_energy_temperature") is not True
        or acceptance.get("SCC_error_detected") is not False
    ):
        raise ValueError("Stage process or production acceptance is not successful")
    if receipt.get("frames_censored") != 0 or receipt.get("frames") != 100:
        raise ValueError("Expected 100 uncensored observables frames")
    roles = receipt.get("roles_1based", {})
    if any(
        roles.get(name) != atom_id
        for name, atom_id in zip(ATOM_NAMES, ATOM_IDS, strict=True)
    ):
        raise ValueError("Observable atom roles differ from sparse-coordinate mapping")

    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError("Observable CSV is missing a header")
        old_rows = list(reader)
        columns = list(reader.fieldnames)
    required = {"time_ps", "q_h", *OLD_DISTANCE_COLUMNS}
    if not required <= set(columns) or any(
        name in columns for name in ("r_nalpha_hg1", "q_act", "q_del")
    ):
        raise ValueError("Observable CSV has missing or already-augmented columns")

    coordinate_path = trajectory
    indices = np.asarray(ATOM_IDS, dtype=np.int64) - 1
    sparse_frames: list[np.ndarray] = []
    times: list[float] = []
    boxes: list[np.ndarray] = []
    with netcdf_file(str(coordinate_path), "r", mmap=True) as nc:
        coordinates = nc.variables["coordinates"]
        time_var = nc.variables["time"]
        lengths = nc.variables["cell_lengths"]
        angles = nc.variables["cell_angles"]
        if coordinates.shape[0] != len(old_rows) or coordinates.shape[1] <= int(
            indices.max()
        ):
            raise ValueError(
                "Native trajectory frame or atom count differs from CSV/mapping"
            )
        if (
            time_var.shape != (len(old_rows),)
            or lengths.shape[0] != len(old_rows)
            or angles.shape[0] != len(old_rows)
        ):
            raise ValueError(
                "Native trajectory time/cell frame counts do not match CSV"
            )
        for frame_index in range(len(old_rows)):
            sparse_frames.append(
                np.array(coordinates[frame_index, indices, :], copy=True)
            )
            times.append(float(time_var[frame_index]))
            boxes.append(
                np.r_[
                    np.array(lengths[frame_index], copy=True),
                    np.array(angles[frame_index], copy=True),
                ]
            )
        del coordinates, time_var, lengths, angles
    augmented = augment_arrays(
        np.stack(sparse_frames), np.asarray(times), np.stack(boxes), old_rows
    )

    output.mkdir(parents=True, exist_ok=False)
    new_columns = ["r_nalpha_hg1", "q_act", "q_del"]
    with (output / "observables.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns + new_columns)
        writer.writeheader()
        for row, extra in zip(old_rows, augmented, strict=True):
            writer.writerow({**row, **extra})
    source_paths = (
        process_path,
        acceptance_path,
        receipt_path,
        csv_path,
        trajectory,
    )
    source_sha256 = {
        str(path.resolve()): raw_hash if path == trajectory else _sha256(path)
        for path in source_paths
    }
    distance_errors = {
        name: max(
            abs(item[name] - float(row[name]))
            for row, item in zip(old_rows, augmented, strict=True)
        )
        for name in OLD_DISTANCE_COLUMNS
    }
    q_h_errors = [
        abs(item["q_act"] + item["q_del"] - float(row["q_h"]))
        for row, item in zip(old_rows, augmented, strict=True)
    ]
    provenance = {
        "stage": str(stage),
        "source_sha256": source_sha256,
        "csv_sha256": _sha256(output / "observables.csv"),
        "rows_count": len(old_rows),
        "source_shape": {
            "frames": len(old_rows),
            "sparse_atoms": 5,
            "coordinates": [len(old_rows), 5, 3],
        },
        "atom_ids_1based": dict(zip(ATOM_NAMES, ATOM_IDS, strict=True)),
        "frames_censored": 0,
        "distance_tolerance_A": 1e-8,
        "time_tolerance_ps": 1e-6,
        "parity": {
            "distance_columns_max_abs_error_A": distance_errors,
            "q_h_max_abs_error_A": max(q_h_errors),
            "all_four_distances_passed": True,
            "q_h_partition_passed": True,
        },
        "pbc_method": "neb_umbrella_observables._distances minimum-image implementation",
    }
    (output / "augmentation.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance
