"""Box-aware per-frame CV, proton-pool and cross-bias diagnostics."""

from __future__ import annotations

import csv
import json
import logging
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask
from scipy.io import netcdf_file

from nyl12_covalent.amber_md_records import instantaneous_md_text
from nyl12_covalent.dftb_scc_validation import audit_dftb_scc
from nyl12_covalent.neb_umbrella_observables import extract_observables
from nyl12_covalent.relay_umbrella_artifacts import digest
from nyl12_covalent.relay_umbrella_cv import (
    RelayBiasWindow,
    amber_bias_energy,
    cross_bias_energies,
    raw_features,
)

LOGGER = logging.getLogger(__name__)
BOLTZMANN_KCAL_MOL_K = 0.00198720425864083
DT_PS = 0.0005
ROLES = {
    "carbon": 19916,
    "leaving_n": 19929,
    "og1": 13916,
    "hg1": 13917,
    "nalpha": 13905,
    "nalpha_h1": 13906,
    "nalpha_h2": 13907,
}
NUMERIC_DIAGNOSTIC_COLUMNS = (
    "r_c_n",
    "r_c_og1",
    "r_og1_hg1",
    "r_nleave_hg1",
    "cn_nleave_hpool",
    "cn_og1_hpool",
    "cn_progress",
    "cn_nalpha_hpool",
    "r_nalpha_h1",
    "r_nalpha_h2",
)
OWNER_CATEGORICAL_COLUMNS = (
    "hg1_owner",
    "nalpha_h1_owner",
    "nalpha_h2_owner",
    "owner_changed",
)
OWNER_STRING_COLUMNS = (
    "changed_hydrogen_ids",
    "hydrogen_owners",
    "hydrogen_owner_distances_A",
    "ambiguous_hydrogens_within_1_30_A",
)


def _read_native_features(
    trajectory: Path, step_origin: int, steps: int, output_stride: int
) -> tuple[list[dict[str, float]], float]:
    """Read all full-coordinate frames with their native periodic cell."""
    records = []
    maximum_coordinate = 0.0
    with netcdf_file(str(trajectory), "r", mmap=False) as stream:
        coordinates = stream.variables["coordinates"]
        times = stream.variables["time"]
        lengths = stream.variables["cell_lengths"]
        angles = stream.variables["cell_angles"]
        if len(coordinates.shape) != 3 or coordinates.shape[2] != 3:
            raise ValueError(
                "Trajectory coordinates must have shape (frames, atoms, 3)"
            )
        if coordinates.shape[0] == 0:
            raise ValueError("Production trajectory has no frames")
        time_values = np.asarray(times[:])
        integration_steps = _validate_trajectory_times(
            time_values, step_origin, steps, output_stride
        )
        for index in range(coordinates.shape[0]):
            xyz = np.asarray(coordinates[index], dtype=float)
            maximum_coordinate = max(maximum_coordinate, float(np.max(np.abs(xyz))))
            box = np.r_[lengths[index], angles[index]]
            record = raw_features(xyz, box)
            record["frame"] = index + 1
            record["time_ps"] = float(times[index])
            record["step"] = integration_steps[index]
            if not all(math.isfinite(v) for v in record.values()):
                raise ValueError(
                    f"Nonfinite CV or time in trajectory frame {index + 1}"
                )
            records.append(record)
    return records, 2.0 * math.sqrt(3.0) * abs(
        float(np.spacing(np.float32(maximum_coordinate)))
    )


def _validate_trajectory_times(
    time_values: np.ndarray,
    step_origin: int,
    steps: int,
    output_stride: int,
) -> list[int]:
    """Validate native trajectory timing against the scheduled output grid.

    The time variable's own floating-point ULP bounds serialization error. The
    sampled step sequence, frame count, and monotonicity remain exact gates.
    """
    if time_values.ndim != 1 or not np.issubdtype(time_values.dtype, np.floating):
        raise ValueError("Trajectory times must be a one-dimensional float array")
    if steps <= 0 or output_stride <= 0 or steps % output_stride:
        raise ValueError("Production steps must be divisible by trajectory stride")
    expected_count = steps // output_stride
    if len(time_values) != expected_count:
        raise ValueError(
            f"Trajectory has {len(time_values)} frames; expected {expected_count}"
        )
    serialized = time_values.astype(np.float64)
    if not np.all(np.isfinite(serialized)):
        raise ValueError("Trajectory times must be finite")
    if np.any(np.diff(serialized) <= 0.0):
        raise ValueError("Trajectory times must be strictly increasing")

    local_steps: list[int] = []
    for index, time_ps in enumerate(serialized):
        expected_local_step = (index + 1) * output_stride
        global_step = time_ps / DT_PS
        local_step = global_step - step_origin
        rounded_step = int(round(local_step))
        native_time = np.asarray(time_values[index], dtype=time_values.dtype)
        half_ulp_steps = abs(float(np.spacing(native_time))) / (2.0 * DT_PS)
        float64_allowance = (
            8.0
            * np.finfo(np.float64).eps
            * max(abs(global_step), abs(step_origin), expected_local_step, 1.0)
        )
        if (
            rounded_step != expected_local_step
            or abs(local_step - expected_local_step)
            > half_ulp_steps + float64_allowance
        ):
            raise ValueError(
                f"Trajectory frame {index + 1} is off the integration grid"
            )
        local_steps.append(rounded_step)
    return local_steps


def _step_origin_from_mdin(mdin: str, restart_time_ps: float) -> int:
    """Derive the native trajectory clock origin from Amber's explicit irest.

    Args:
        mdin: Contents of the stage's Amber input file.
        restart_time_ps: Time stored in the stage input restart file.

    Returns:
        Integration-step origin for native NetCDF times.

    Raises:
        ValueError: If irest is missing, duplicated, unsupported, or invalid.
    """
    uncommented = "\n".join(line.split("!", 1)[0] for line in mdin.splitlines())
    declarations = re.findall(r"\birest\s*=\s*([^,\s/]+)", uncommented, re.IGNORECASE)
    if len(declarations) != 1:
        raise ValueError("Amber input must declare irest exactly once")
    try:
        irest = int(declarations[0])
    except ValueError as error:
        raise ValueError("Amber input irest must be integer 0 or 1") from error
    if irest == 0:
        return 0
    if irest != 1:
        raise ValueError(f"Unsupported Amber irest value: {irest}")
    if not math.isfinite(restart_time_ps) or restart_time_ps < 0.0:
        raise ValueError("Restart time must be finite and nonnegative")
    return int(round(restart_time_ps / DT_PS))


def _read_dumpave(path: Path) -> dict[int, tuple[float, float]]:
    """Read native Amber step and two restraint CVs, preserving native precision."""
    values: dict[int, tuple[float, float]] = {}
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 3:
            raise ValueError(
                f"Expected step and two DUMPAVE CVs at {path}:{line_number}"
            )
        try:
            step = int(fields[0])
            cv_values = (float(fields[1]), float(fields[2]))
        except ValueError as error:
            raise ValueError(
                f"Invalid DUMPAVE record at {path}:{line_number}"
            ) from error
        if step == 0 and values.get(0) == cv_values:
            LOGGER.info(
                "Coalesced identical unsampled startup CV record at %s:%d",
                path,
                line_number,
            )
            continue
        if step in values or not all(math.isfinite(value) for value in cv_values):
            raise ValueError(f"Duplicate/nonfinite DUMPAVE record at step {step}")
        values[step] = cv_values
    if not values:
        raise ValueError(f"No native CV records in {path}")
    return values


def _native_bias_records(
    feature_records: list[dict[str, float]],
    dumpave: Mapping[int, tuple[float, float]],
    window: RelayBiasWindow,
    serialization_tolerance_A: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Evaluate bias from native DUMPAVE CVs and check them against NetCDF."""
    count = len(feature_records)
    bias = np.full(count, np.nan)
    native_cvs = np.full((count, 2), np.nan)
    bounds = np.full(count, np.nan)
    evaluated = np.zeros(count, dtype=bool)
    epsilon_A = 0.0005 + serialization_tolerance_A
    for index, record in enumerate(feature_records):
        step = int(record["step"])
        if step not in dumpave:
            if index != count - 1:
                raise ValueError(
                    f"Native DUMPAVE lacks nonterminal trajectory step {step}"
                )
            continue
        values = dumpave[step]
        native_cvs[index] = values
        for name, native_value in zip(window.feature_names, values, strict=True):
            if abs(record[name] - native_value) > epsilon_A:
                raise ValueError(
                    f"Native DUMPAVE/NetCDF CV mismatch at step {step}: {name}"
                )
        bias[index] = sum(
            amber_bias_energy(value, center, force, window.outer_extent_A)
            for value, center, force in zip(
                values, window.centers_A, window.force_constants, strict=True
            )
        )
        bounds[index] = (
            sum(
                2.0
                * force
                * min(abs(value - center) + epsilon_A, window.outer_extent_A)
                * epsilon_A
                for value, center, force in zip(
                    values, window.centers_A, window.force_constants, strict=True
                )
            )
            + 1e-4
        )
        evaluated[index] = True
    if evaluated.sum() != count - 1 or evaluated[-1]:
        raise ValueError(
            "Only the final trajectory frame may lack a native DUMPAVE record"
        )
    return bias, native_cvs, bounds, evaluated, epsilon_A


def _native_stage_receipt(
    stage: Path,
    steps: int,
    temperature_K: float,
    temperature_gate_K: tuple[float, float],
) -> dict[str, Any]:
    """Confirm process completion and DFTB SCC resolution before diagnostics."""
    process = json.loads((stage / "process.json").read_text())
    if process.get("exit_status") != 0 or process.get("ended_utc") is None:
        raise ValueError("Diagnostics require a closed successful native process")
    text = (stage / "run.mdout").read_text()
    scc = audit_dftb_scc(text)
    if scc.get("passed") is not True or scc.get("verified_scc_history") is not True:
        raise ValueError("Diagnostics require resolved DFTB SCC history")
    rows = re.findall(
        r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([-+\d.EeDd]+)\s+TEMP\(K\)\s*=\s*([-+\d.EeDd]+)",
        instantaneous_md_text(text),
    )
    if not rows or int(rows[-1][0]) != steps:
        raise ValueError("Native stage did not complete its declared step count")
    final_temperature = float(rows[-1][2].replace("D", "E"))
    if (
        not math.isfinite(final_temperature)
        or not temperature_gate_K[0] <= final_temperature <= temperature_gate_K[1]
    ):
        raise ValueError(
            "Final production temperature is outside the declared diagnostic gate"
        )
    return {
        "exit_status": process["exit_status"],
        "electronic_history": scc,
        "final_temperature_K": final_temperature,
    }


def _write_features_csv(
    output: Path,
    observables_csv: Path,
    features: list[dict[str, float]],
    native_bias: np.ndarray,
    native_cvs: np.ndarray,
    native_bias_error_bound: np.ndarray,
    native_bias_evaluated: np.ndarray,
    computed_bias: np.ndarray,
) -> str:
    """Merge all-H observables with three periodic relay CVs and bias values."""
    with observables_csv.open(newline="") as source:
        source_rows = list(csv.DictReader(source))
    if any(
        len(value) != len(features)
        for value in (
            source_rows,
            native_bias,
            native_cvs,
            native_bias_error_bound,
            native_bias_evaluated,
            computed_bias,
        )
    ):
        raise ValueError("Observable, CV and native-bias frame counts differ")
    with output.open("x", newline="") as stream:
        fieldnames = [
            *source_rows[0].keys(),
            "q_acyl_A",
            "q_act_A",
            "q_del_A",
            "q_h_A",
            "native_bias_kcal_mol",
            "native_bias_evaluated",
            "native_cv_evaluated",
            "native_cv_1_A",
            "native_cv_2_A",
            "native_bias_error_bound_kcal_mol",
            "computed_bias_kcal_mol",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for index, (source, feature, bias) in enumerate(
            zip(source_rows, features, native_bias, strict=True)
        ):
            if "owner_changed" in source:
                flag = source["owner_changed"].strip().lower()
                if flag not in {"true", "false", "1", "0"}:
                    raise ValueError(
                        "Invalid proton-owner change flag in observable CSV"
                    )
                source["owner_changed"] = int(flag in {"true", "1"})
            if abs(float(source["q_acyl"]) - feature["q_acyl"]) > 1e-6:
                raise ValueError(
                    "q_acyl differs from the established periodic observable"
                )
            if abs(float(source["q_h"]) - feature["q_h"]) > 1e-6:
                raise ValueError(
                    "Legacy q_H differs from periodic relay feature identity"
                )
            if abs(feature["q_h"] - feature["q_act"] - feature["q_del"]) > 1e-8:
                raise ValueError("q_H != q_act + q_del for this frame")
            writer.writerow(
                {
                    **source,
                    "q_acyl_A": feature["q_acyl"],
                    "q_act_A": feature["q_act"],
                    "q_del_A": feature["q_del"],
                    "q_h_A": feature["q_h"],
                    "native_bias_kcal_mol": float(bias)
                    if native_bias_evaluated[index]
                    else "",
                    "native_bias_evaluated": int(native_bias_evaluated[index]),
                    "native_cv_evaluated": int(native_bias_evaluated[index]),
                    "native_cv_1_A": float(native_cvs[index, 0])
                    if native_bias_evaluated[index]
                    else "",
                    "native_cv_2_A": float(native_cvs[index, 1])
                    if native_bias_evaluated[index]
                    else "",
                    "native_bias_error_bound_kcal_mol": float(
                        native_bias_error_bound[index]
                    )
                    if native_bias_evaluated[index]
                    else "",
                    "computed_bias_kcal_mol": float(computed_bias[index]),
                }
            )
    return digest(output)


def block_diagnostics(
    cfg: Mapping[str, Any],
    row: Mapping[str, Any],
    stage: Path,
    topology: Any,
) -> dict[str, Any]:
    """Write periodic raw features, exact all-window bias matrix and proof.

    The execution runner calls this after process completion, SCC acceptance,
    and frozen-shell checkpoint checks. Every valid trajectory frame is kept.
    """
    campaign = stage.resolve().parents[1]
    trajectory = stage / "trajectory.nc"
    steps = next(
        spec["steps"]
        for spec in row["stages"]
        if (campaign / spec["input"]).parent == stage.resolve()
    )
    native_acceptance = json.loads((stage / "acceptance.json").read_text())
    if (
        native_acceptance.get("completed_steps") != steps
        or native_acceptance.get("every_integration_step_SCC_converged") is not True
        or native_acceptance.get("finite_energy_temperature") is not True
        or not native_acceptance.get("checkpoint_sha256")
    ):
        raise ValueError("Diagnostics require immutable native acceptance")
    if (campaign / "prepared.sha256").read_text().strip() != digest(
        campaign / "prepared.json"
    ):
        raise ValueError("Prepared manifest hash sidecar mismatch")
    gate = tuple(cfg.get("diagnostic_temperature_gate_K", (225.0, 375.0)))
    if len(gate) != 2 or gate[0] >= gate[1]:
        raise ValueError("Invalid diagnostic temperature gate")
    if cfg.get("production_data", True) is not True and gate != (225.0, 375.0):
        declared = cfg.get("diagnostic_temperature_gate_K")
        if declared is None:
            raise ValueError("A nonproduction fixture gate must be explicit")
    elif cfg.get("production_data", True) is True and gate != (225.0, 375.0):
        raise ValueError("Production diagnostics require the fixed 300 K +/-75 K gate")
    native = _native_stage_receipt(
        stage, steps, float(cfg["temperature_K"]), (float(gate[0]), float(gate[1]))
    )
    mdin = (stage / "run.mdin").read_text()
    input_restart = Rst7.open(str(stage / "input.ncrst"))
    step_origin = _step_origin_from_mdin(mdin, float(input_restart.time))
    stride_match = re.search(r"\bntwx\s*=\s*(\d+)", mdin, re.IGNORECASE)
    if stride_match is None:
        raise ValueError("Production input does not declare trajectory stride ntwx")
    output_stride = int(stride_match.group(1))
    feature_records, serialization_tolerance_A = _read_native_features(
        trajectory, step_origin, steps, output_stride
    )
    source_manifest = json.loads((campaign / "prepared.json").read_text())
    qm_mask = source_manifest["system"]["qm_atom_ids_mask"]
    qm_ids = [i + 1 for i in AmberMask(topology, qm_mask).Selected()]
    observables_dir = stage / "observables"
    extract_observables(
        campaign / cfg["topology"],
        trajectory,
        ROLES,
        qm_ids,
        observables_dir,
    )
    with (observables_dir / "observables.csv").open(newline="") as stream:
        observable_rows = list(csv.DictReader(stream))
    if len(observable_rows) != len(feature_records):
        raise ValueError("Full-pool and raw-CV frame counts differ")

    windows = [RelayBiasWindow(**candidate["bias"]) for candidate in cfg["windows"]]
    feature_arrays = {
        name: np.asarray([record[name] for record in feature_records], dtype=float)
        for name in ("q_acyl", "q_act", "q_del", "q_h")
    }
    cross_bias = cross_bias_energies(feature_arrays, windows)
    window_index = next(
        i
        for i, candidate in enumerate(cfg["windows"])
        if candidate["label"] == row["label"]
    )
    sampled_window = windows[window_index]
    dumpave = _read_dumpave(stage / "umbrella_coordinate.dat")
    native_bias, native_cvs, native_bound, native_evaluated, epsilon_A = (
        _native_bias_records(
            feature_records, dumpave, sampled_window, serialization_tolerance_A
        )
    )
    bias_for_sampled = cross_bias[window_index]
    csv_path = stage / "raw_features.csv"
    csv_sha = _write_features_csv(
        csv_path,
        observables_dir / "observables.csv",
        feature_records,
        native_bias,
        native_cvs,
        native_bound,
        native_evaluated,
        bias_for_sampled,
    )

    beta = 1.0 / (BOLTZMANN_KCAL_MOL_K * float(cfg["temperature_K"]))
    npz_path = stage / "samples_and_cross_bias.npz"
    np.savez_compressed(
        npz_path,
        feature_names=np.asarray(["q_acyl", "q_act", "q_del"]),
        feature_values=np.column_stack(
            [feature_arrays[name] for name in ("q_acyl", "q_act", "q_del")]
        ),
        q_h_A=feature_arrays["q_h"],
        window_labels=np.asarray([candidate["label"] for candidate in cfg["windows"]]),
        window_feature_names=np.asarray(
            [candidate["biased_cv_names"] for candidate in cfg["windows"]]
        ),
        bias_kcal_mol=cross_bias,
        u_kn_bias_only=beta * cross_bias,
    )
    process_path = stage / "process.json"
    process_sha = digest(process_path)
    native_artifacts = {
        str(path.relative_to(campaign)): digest(path)
        for path in (
            trajectory,
            stage / "restart.ncrst",
            stage / "run.mdout",
            stage / "umbrella_coordinate.dat",
            observables_dir / "observables.csv",
            observables_dir / "observables.json",
        )
    }
    input_identity = json.loads((stage / "input_identity.json").read_text())
    sourceproof = {
        "schema_version": 1,
        "window": row["label"],
        "prepared_sha256": digest(campaign / "prepared.json"),
        "prepared_sha256_sidecar": (campaign / "prepared.sha256").read_text().strip(),
        "analysis_request_sha256": digest(campaign / "analysis_request.json"),
        "physical_hamiltonian_sha256": cfg["physical_hamiltonian_sha256"],
        "source_evidence_sha256": {
            name: value
            for name, value in cfg["frozen_files_sha256"].items()
            if name.startswith("source_evidence/")
        },
        "source_snapshot_sha256": {
            name: value
            for name, value in cfg["frozen_files_sha256"].items()
            if name.startswith("source_snapshot/")
        },
        "input_identity": input_identity,
        "acceptance": native_acceptance,
        "native_artifacts_sha256": native_artifacts,
        "native_bias_source": "amber_dumpave_independent_CVs",
        "native_bias_feature_names": list(sampled_window.feature_names),
        "native_cv_quantization_A": 0.0005,
        "coordinate_serialization_tolerance_A": serialization_tolerance_A,
        "native_cv_total_tolerance_A": epsilon_A,
        "native_bias_evaluated_frames": int(native_evaluated.sum()),
        "native_bias_unevaluated_terminal_frames": int((~native_evaluated).sum()),
    }
    sourceproof_path = stage / "sourceproof.json"
    sourceproof_path.write_text(
        json.dumps(sourceproof, indent=2, sort_keys=True) + "\n"
    )
    source_sha = {
        "sourceproof.json": digest(sourceproof_path),
        "umbrella_coordinate.dat": digest(stage / "umbrella_coordinate.dat"),
    }
    duration_ps = steps * DT_PS
    proof = {
        "accepted": True,
        "exit_status": native["exit_status"],
        "frames_censored": 0,
        "temperature_K": float(cfg["temperature_K"]),
        "final_temperature_K": native["final_temperature_K"],
        "physical_hamiltonian_sha256": cfg["physical_hamiltonian_sha256"],
        "csv_sha256": csv_sha,
        "process_sha256": process_sha,
        "source_sha256": source_sha,
        "native_artifacts_sha256": native_artifacts,
        "native_bias_source": "amber_dumpave_independent_CVs",
        "native_cv_quantization_A": 0.0005,
        "coordinate_serialization_tolerance_A": serialization_tolerance_A,
        "native_cv_total_tolerance_A": epsilon_A,
        "native_bias_evaluated_frames": int(native_evaluated.sum()),
        "native_bias_unevaluated_terminal_frames": int((~native_evaluated).sum()),
        "window": row["label"],
        "frames": len(feature_records),
        "sampled_duration_ps": duration_ps,
        "csv_file": "raw_features.csv",
        "process_file": "process.json",
        "coordinates_finite": True,
        "energies_finite": True,
        "electronic_history": native["electronic_history"],
    }
    proof_path = stage / "analysis_acceptance.json"
    proof_path.write_text(json.dumps(proof, indent=2, sort_keys=True) + "\n")
    analysis_request = {
        "windows": cfg["windows"],
        "path_3cv_A": cfg["path_3cv_A"],
        "temperature_K": cfg["temperature_K"],
        "physical_hamiltonian_sha256": cfg["physical_hamiltonian_sha256"],
        "accepted_block_receipts": [str(proof_path.relative_to(campaign))],
        "es_search_bounds": [-0.10, 0.10],
        "acyl_search_bounds": [0.90, 1.10],
        "diagnostic_columns": list(NUMERIC_DIAGNOSTIC_COLUMNS),
        "owner_categorical_columns": list(OWNER_CATEGORICAL_COLUMNS),
        "owner_string_columns": list(OWNER_STRING_COLUMNS),
        "frames_retained": len(feature_records),
        "frames_censored": 0,
        "bias_matrix_shape": list(cross_bias.shape),
        "csv_sha256": csv_sha,
        "cross_bias_npz_sha256": digest(npz_path),
        "bias_convention": "Amber k*(xi-center)^2; no half factor; linear tails outside center +/- 10 A.",
    }
    (stage / "analysis_request.json").write_text(
        json.dumps(analysis_request, indent=2, sort_keys=True) + "\n"
    )
    return proof
