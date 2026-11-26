"""Validate immutable direct-US continuation plans and requested native stages."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

from .relay_umbrella_analysis_io import digest

BLOCK_PS = 5.0
FIRST_EXTENSION_BLOCK = 7
LAST_EXTENSION_BLOCK = 12
_TRAJECTORY_HASH_CACHE: dict[Path, tuple[int, int, int, int, int, str]] = {}


def _trajectory_digest(path: Path) -> str:
    """Hash a native trajectory once per unchanged filesystem identity."""
    resolved = path.resolve(strict=True)
    stat = resolved.stat()
    identity = (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )
    cached = _TRAJECTORY_HASH_CACHE.get(resolved)
    if cached is not None and cached[:5] == identity:
        return cached[5]
    value = digest(resolved)
    _TRAJECTORY_HASH_CACHE[resolved] = (*identity, value)
    return value


def requested_block_count(horizon_ps: float) -> int:
    """Return required absolute production block count for a supported horizon.

    Args:
        horizon_ps: Total native production duration per window, in ps.

    Returns:
        Number of 5 ps production blocks needed through ``horizon_ps``.

    Raises:
        ValueError: If the horizon is not an extension multiple from 35 to 60 ps.
    """
    if (
        not math.isfinite(horizon_ps)
        or horizon_ps < 35.0
        or horizon_ps > 60.0
        or not math.isclose(horizon_ps / BLOCK_PS, round(horizon_ps / BLOCK_PS))
    ):
        raise ValueError("Continuation horizon must be 35, 40, 45, 50, 55, or 60 ps")
    return int(round(horizon_ps / BLOCK_PS))


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def validate_continuation_plan(
    campaign: Path,
    continuation_plan: Path,
    horizon_ps: float,
    augmented_root: Path | None = None,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, str]]:
    """Validate the plan and materialized stage inputs through one horizon.

    Future stage acceptance is deliberately outside this validator: only native
    input identities required by the requested horizon are checked here.

    Args:
        campaign: Original direct-US root containing ``prepared.json``.
        continuation_plan: Immutable ``plan.json`` produced by the native prep.
        horizon_ps: Requested total production horizon, from 35 through 60 ps.
        augmented_root: Optional derived CSV root to verify against the plan.

    Returns:
        Per-window requested stage records and hashes for provenance.

    Raises:
        FileNotFoundError: If a required plan or requested stage input is absent.
        ValueError: If plan, prepared, native stage, or hash identities disagree.
    """
    requested_blocks = requested_block_count(horizon_ps)
    campaign = Path(campaign).resolve(strict=True)
    plan_path = Path(continuation_plan).resolve(strict=True)
    plan_hash_path = plan_path.with_name("plan.sha256")
    expected_plan_hash = plan_hash_path.read_text().strip()
    if digest(plan_path) != expected_plan_hash:
        raise ValueError("Continuation plan digest mismatch")
    plan = _read_json(plan_path)
    prepared_path = campaign / "prepared.json"
    if Path(plan.get("campaign_root", "")).resolve() != campaign:
        raise ValueError("Continuation plan campaign root mismatch")
    if Path(plan.get("prepared", "")).resolve() != prepared_path.resolve():
        raise ValueError("Continuation plan does not reference this prepared.json")
    if (
        augmented_root is not None
        and Path(plan.get("augmentation_root", "")).resolve()
        != Path(augmented_root).resolve()
    ):
        raise ValueError("Continuation augmentation root differs from requested input")
    if digest(prepared_path) != plan.get("prepared_sha256"):
        raise ValueError("Prepared campaign differs from continuation plan")
    prepared = _read_json(prepared_path)
    expected_windows = [row["label"] for row in prepared["windows"]]
    labels = [row.get("label") for row in plan.get("windows", [])]
    if len(expected_windows) != 63 or labels != expected_windows:
        raise ValueError("Continuation plan window order differs from campaign")

    source_hashes = plan.get("source_sha256")
    if not isinstance(source_hashes, dict):
        raise ValueError("Continuation plan lacks source_sha256 identities")
    sources = {
        str(plan_path): digest(plan_path),
        str(plan_hash_path): digest(plan_hash_path),
        str(prepared_path): digest(prepared_path),
    }
    for filename, expected in source_hashes.items():
        path = Path(filename)
        actual = digest(path)
        if actual != expected:
            raise ValueError(f"Continuation source hash mismatch: {path}")
        sources[str(path)] = actual

    flat = plan.get("stages")
    if not isinstance(flat, list):
        raise ValueError("Continuation plan lacks a flat stage manifest")
    per_window = {label: [] for label in expected_windows}
    seen: set[tuple[str, int]] = set()
    stage_by_window = {row["label"]: row.get("stages", []) for row in plan["windows"]}
    for record in flat:
        label, block = record.get("window"), record.get("block")
        if label not in per_window or type(block) is not int:
            raise ValueError("Invalid continuation window/block identity")
        if not FIRST_EXTENSION_BLOCK <= block <= LAST_EXTENSION_BLOCK:
            raise ValueError("Continuation block must be 07 through 12")
        pair = (label, block)
        if pair in seen:
            raise ValueError(f"Duplicate continuation stage: {pair}")
        seen.add(pair)
        expected_label = f"03_production_block_{block:02d}_5ps"
        if record.get("label") != expected_label:
            raise ValueError(f"Continuation label/block mismatch: {record}")
        window_stages = stage_by_window[label]
        if not any(
            item.get("label") == expected_label
            and item.get("kind") == "production"
            and item.get("block") == block
            for item in window_stages
        ):
            raise ValueError(f"Per-window stage list omits {label}/{expected_label}")
        stage_path = Path(record.get("stage", "")).resolve()
        expected_stage = (campaign / label / expected_label).resolve()
        if stage_path != expected_stage:
            raise ValueError(f"Continuation stage path mismatch: {stage_path}")
        if block <= requested_blocks:
            input_hashes = {
                "run.mdin": record.get("mdin_sha256"),
                "umbrella.RST": record.get("restraint_sha256"),
            }
            for name, expected in input_hashes.items():
                path = stage_path / name
                if not expected or digest(path) != expected:
                    raise ValueError(f"Continuation stage input hash mismatch: {path}")
                sources[str(path)] = expected
            per_window[label].append(record)

    expected_pairs = {
        (label, block)
        for label in expected_windows
        for block in range(FIRST_EXTENSION_BLOCK, requested_blocks + 1)
    }
    if not expected_pairs <= seen:
        missing = sorted(expected_pairs - seen)
        raise ValueError(
            f"Continuation plan skips required production blocks: {missing[:5]}"
        )
    for label in per_window:
        per_window[label].sort(key=lambda record: record["block"])
        if len(per_window[label]) != requested_blocks - FIRST_EXTENSION_BLOCK + 1:
            raise ValueError(f"Incomplete requested continuation stages for {label}")
    return per_window, sources


def validate_accepted_stage(
    record: dict[str, Any], augmented_root: Path
) -> dict[str, str]:
    """Require native acceptance and coordinate augmentation for one stage.

    Args:
        record: Continuation stage record from the verified flat manifest.
        augmented_root: Derived augmentation root named by the plan.

    Returns:
        Hashes of validated native compact files and derived outputs.

    Raises:
        FileNotFoundError: If the accepted native or derived result is absent.
        ValueError: If acceptance, frame coverage, or source identities fail.
    """
    stage = Path(record["stage"])
    native_paths = {
        "process": stage / "process.json",
        "acceptance": stage / "acceptance.json",
        "observables": stage / "observables/observables.json",
        "native_csv": stage / "observables/observables.csv",
        "trajectory": stage / "trajectory.nc",
    }
    hashes = {
        key: _trajectory_digest(path) if key == "trajectory" else digest(path)
        for key, path in native_paths.items()
    }
    process = _read_json(native_paths["process"])
    acceptance = _read_json(native_paths["acceptance"])
    native_receipt = _read_json(native_paths["observables"])
    if process.get("exit_status") != 0 or not process.get("ended_utc"):
        raise ValueError(f"Continuation process is not accepted: {stage}")
    if (
        acceptance.get("completed_steps") != 10000
        or acceptance.get("finite_energy_temperature") is not True
        or acceptance.get("SCC_error_detected") is not False
    ):
        raise ValueError(f"Continuation acceptance failed: {stage}")
    if (
        native_receipt.get("frames") != 100
        or native_receipt.get("frames_censored") != 0
        or native_receipt.get("csv_sha256") != digest(native_paths["native_csv"])
    ):
        raise ValueError(f"Continuation native observable receipt failed: {stage}")
    trajectory_expected = native_receipt.get("source_sha256", {}).get(
        str(native_paths["trajectory"])
    )
    if not trajectory_expected or trajectory_expected != _trajectory_digest(
        native_paths["trajectory"]
    ):
        raise ValueError(f"Continuation raw trajectory hash failed: {stage}")
    with native_paths["native_csv"].open(newline="") as stream:
        if sum(1 for _ in csv.DictReader(stream)) != 100:
            raise ValueError(f"Continuation native CSV frame count failed: {stage}")

    derived = Path(augmented_root) / record["window"] / record["label"]
    derived_receipt_path = derived / "augmentation.json"
    derived_csv = derived / "observables.csv"
    derived_receipt = _read_json(derived_receipt_path)
    if (
        derived_receipt.get("rows_count") != 100
        or derived_receipt.get("frames_censored") != 0
        or derived_receipt.get("csv_sha256") != digest(derived_csv)
        or derived_receipt.get("parity", {}).get("all_four_distances_passed")
        is not True
        or derived_receipt.get("parity", {}).get("q_h_partition_passed") is not True
    ):
        raise ValueError(f"Continuation augmentation receipt failed: {derived}")
    with derived_csv.open(newline="") as stream:
        if sum(1 for _ in csv.DictReader(stream)) != 100:
            raise ValueError(
                f"Continuation augmented CSV frame count failed: {derived}"
            )
    source_hashes = derived_receipt.get("source_sha256", {})
    for key, path in (
        ("process.json", native_paths["process"]),
        ("acceptance.json", native_paths["acceptance"]),
        ("observables/observables.json", native_paths["observables"]),
        ("observables/observables.csv", native_paths["native_csv"]),
        ("trajectory.nc", native_paths["trajectory"]),
    ):
        expected = source_hashes.get(str(path))
        if not expected or expected != digest(path):
            raise ValueError(f"Continuation augmentation source mismatch: {path}")
        hashes[str(path)] = (
            _trajectory_digest(path) if key == "trajectory.nc" else expected
        )
    hashes[str(derived_receipt_path)] = digest(derived_receipt_path)
    hashes[str(derived_csv)] = digest(derived_csv)
    return hashes
