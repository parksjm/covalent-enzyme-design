"""Manifest, reuse, and artifact validation for local design campaigns."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    """Hash a file without changing it."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    """Atomically publish a machine-readable status or provenance record."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def case_path(campaign: Path, job: dict[str, Any]) -> Path:
    """Return the new model directory using the stable one-based seed folder."""
    return (
        campaign / "models" / job["variant"] / f"seed_{int(job['model_index']) + 1:02d}"
    )


def component_path(campaign: Path, job: dict[str, Any]) -> Path:
    """Return the per-prediction Boltz component graph pickle."""
    variant, index = job["variant"], int(job["model_index"])
    return (
        campaign
        / "predictions/runs"
        / variant
        / f"model_{index}"
        / "output"
        / f"boltz_results_{variant}/processed/mols/{variant}.pkl"
    )


def validate_jobs(
    campaign: Path, manifest: dict[str, Any], root: Path
) -> list[dict[str, Any]]:
    """Validate new-job lineage and exact agreement with prepared/jobs.csv."""
    jobs, reused = manifest.get("jobs"), manifest.get("reused_models")
    if not isinstance(jobs, list) or not 1 <= len(jobs) <= 52:
        raise ValueError("candidate_manifest.jobs must contain 1–52 new models")
    if not isinstance(reused, list):
        raise ValueError("candidate_manifest.reused_models must be a list")
    keys: set[tuple[str, int]] = set()
    for job in jobs:
        required = {
            "variant",
            "model_index",
            "seed",
            "yaml_path",
            "expected_alpha_sequence",
            "expected_beta_sequence",
        }
        if not required <= job.keys():
            raise ValueError(f"New job lacks fields: {sorted(required - job.keys())}")
        key = (job["variant"], int(job["model_index"]))
        if key in keys:
            raise ValueError(f"Duplicate new model lineage: {key}")
        keys.add(key)
        if (
            len(job["expected_alpha_sequence"]) != 232
            or len(job["expected_beta_sequence"]) != 96
        ):
            raise ValueError(f"Unexpected protein sequence lengths for {key}")
        path = Path(job["yaml_path"])
        path = path if path.is_absolute() else (root / path).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
    new_wt = {int(job["model_index"]) for job in jobs if job["variant"] == "WT"}
    reused_wt = {
        int(item["model_index"]) for item in reused if item.get("variant") == "WT"
    }
    if new_wt != {1, 3} or reused_wt != {0, 2, 4} or len(reused) != 3:
        raise ValueError("Expected new WT indices 1/3 and reused WT indices 0/2/4")
    csv_path = campaign / "prepared/jobs.csv"
    with csv_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required_columns = {
            "campaign",
            "sequence_name",
            "model_index",
            "seed",
            "yaml_path",
        }
        if not required_columns <= set(reader.fieldnames or ()):
            raise ValueError(
                f"jobs.csv lacks columns: {sorted(required_columns - set(reader.fieldnames or ()))}"
            )
        rows = list(reader)
    if len(rows) != len(jobs):
        raise ValueError("jobs.csv row count differs from new candidate manifest")
    row_keys: set[tuple[str, int]] = set()
    jobs_by_key = {(job["variant"], int(job["model_index"])): job for job in jobs}
    for row in rows:
        key = (row["sequence_name"], int(row["model_index"]))
        if key in row_keys:
            raise ValueError(f"Duplicate jobs.csv lineage: {key}")
        row_keys.add(key)
        job = jobs_by_key.get(key)
        if job is None or int(row["seed"]) != int(job["seed"]):
            raise ValueError(f"jobs.csv lineage mismatch for {key}")
        row_yaml, job_yaml = Path(row["yaml_path"]), Path(job["yaml_path"])
        row_yaml = row_yaml if row_yaml.is_absolute() else (root / row_yaml).resolve()
        job_yaml = job_yaml if job_yaml.is_absolute() else (root / job_yaml).resolve()
        if row_yaml.resolve() != job_yaml.resolve():
            raise ValueError(f"jobs.csv YAML mismatch for {key}")
    if row_keys != keys:
        raise ValueError("jobs.csv model set differs from candidate manifest")
    return jobs


def validate_reused(item: dict[str, Any]) -> dict[str, Any]:
    """Verify one immutable reused model’s raw, relaxed, and acceptance lineage."""
    required = {
        "variant",
        "model_index",
        "seed",
        "raw_path",
        "relaxed_path",
        "acceptance_json",
        "source_sha256",
    }
    if not required <= item.keys():
        raise ValueError(
            f"Reused WT model lacks fields: {sorted(required - item.keys())}"
        )
    if item["variant"] != "WT":
        raise ValueError("Only WT models may be reused in this comparison")
    raw, relaxed, acceptance_path = (
        Path(item[name]) for name in ("raw_path", "relaxed_path", "acceptance_json")
    )
    if any(
        not path.is_absolute() or not path.is_file()
        for path in (raw, relaxed, acceptance_path)
    ):
        raise FileNotFoundError("Reused model paths must be absolute existing files")
    raw_hash = sha256(raw)
    if raw_hash != item["source_sha256"]:
        raise ValueError(f"Reused raw source SHA256 mismatch: {raw}")
    acceptance = json.loads(acceptance_path.read_text())
    if Path(acceptance["endpoint"]).resolve() != relaxed.resolve():
        raise ValueError(f"Reused endpoint path mismatch: {acceptance_path}")
    source_copy = Path(acceptance["source"])
    if not source_copy.is_file() or sha256(source_copy) != raw_hash:
        raise ValueError(
            f"Reused relaxation source does not match raw model: {source_copy}"
        )
    process_path = acceptance_path.with_name("process_acceptance.json")
    process = json.loads(process_path.read_text())
    if process.get("status") != "completed" or acceptance.get("process_exit_code") != 0:
        raise ValueError(f"Reused relaxation process did not complete: {process_path}")
    if not isinstance(acceptance.get("force_converged"), bool):
        raise ValueError(
            f"Missing separate force-convergence result: {acceptance_path}"
        )
    return {
        "acceptance": acceptance,
        "process": process,
        "raw": raw,
        "relaxed": relaxed,
        "acceptance_path": acceptance_path,
        "raw_sha256": raw_hash,
    }


def verify_runtime(
    previous: Path, worker_script: Path, cache: Path, boltz_python: Path
) -> dict[str, str]:
    """Verify runtime records from the completed native 45-model campaign."""
    runtime_path = previous / "dispatch/runtime_versions.json"
    checked_path = previous / "dispatch/checked_runtime_manifest.json"
    plan_path = previous / "dispatch/plan.json"
    coordinator_path = previous / "jobs/coordinator_checked/status.json"
    receipt_path = previous / "verification/checkpoint_receipt.json"
    campaign_status_path = previous / "status.json"
    manifest_path = previous / "prepared/candidate_manifest.json"
    config_path = previous / "configuration/relaxation.json"
    required_files = (
        runtime_path,
        checked_path,
        plan_path,
        coordinator_path,
        receipt_path,
        campaign_status_path,
        manifest_path,
        config_path,
    )
    if any(not path.is_file() for path in required_files):
        raise FileNotFoundError("Completed native campaign provenance is incomplete")
    runtime = json.loads(runtime_path.read_text())
    checked = json.loads(checked_path.read_text())
    plan = json.loads(plan_path.read_text())
    coordinator = json.loads(coordinator_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    campaign_status = json.loads(campaign_status_path.read_text())
    if (
        campaign_status.get("status") != "completed"
        or campaign_status.get("planned_models") != 45
        or campaign_status.get("normalized_models") != 45
        or campaign_status.get("relaxed_models") != 45
        or coordinator.get("exit_status") != 0
    ):
        raise ValueError("Previous native 45-model campaign did not complete cleanly")
    if (
        runtime.get("same_boltz_source") is not True
        or runtime.get("boltz_distribution") != "2.2.1"
        or runtime.get("violet", {}).get("torch") != "2.9.1+cu128"
    ):
        raise ValueError(
            "Previous campaign runtime does not match qualified Violet build"
        )
    input_hashes = coordinator.get("input_sha256", {})
    bound_inputs = (checked_path, manifest_path, config_path)
    for path in bound_inputs:
        expected = input_hashes.get(str(path.resolve()))
        if expected is None or sha256(path) != expected:
            raise ValueError(f"Native coordinator input hash mismatch: {path}")
    command = coordinator.get("command", [])
    if "--campaign" not in command or str(previous.resolve()) not in command:
        raise ValueError("Coordinator status is not bound to the previous campaign")
    snapshot = (previous / "source_snapshot").resolve()
    verified_hashes = {str(path.resolve()): sha256(path) for path in required_files}
    for relative_path, expected in checked["sha256"].items():
        path = (previous / relative_path).resolve()
        if not path.is_file() or not path.resolve().is_relative_to(snapshot):
            raise ValueError(
                f"Runtime source is missing or outside frozen snapshot: {path}"
            )
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"Runtime source hash mismatch: {path}")
        verified_hashes[str(path.resolve())] = actual
    for relative_path, expected in plan.get("source_sha256", {}).items():
        path = Path(relative_path)
        if not path.is_file() or not path.resolve().is_relative_to(snapshot):
            raise ValueError(f"Planned source is missing or outside snapshot: {path}")
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"Planned source hash mismatch: {path}")
        verified_hashes[str(path.resolve())] = actual
    if not (snapshot / "run_boltz_cli.py").is_file():
        raise FileNotFoundError(snapshot / "run_boltz_cli.py")
    if not worker_script.is_file() or not cache.is_dir() or not boltz_python.is_file():
        raise FileNotFoundError("Boltz worker, cache, or Python executable is missing")
    checkpoint = cache / "boltz2_conf.ckpt"
    checkpoint_hash = sha256(checkpoint)
    runtime_hash = sha256(runtime_path)
    recorded_runtime = receipt["prediction_runtime"]
    if (
        checkpoint_hash != runtime.get("checkpoint_sha256")
        or checkpoint_hash != recorded_runtime.get("checkpoint_sha256")
        or checkpoint.stat().st_size != runtime.get("checkpoint_bytes")
        or recorded_runtime.get("source_sha256") != runtime_hash
    ):
        raise ValueError("Boltz checkpoint/runtime receipt mismatch")
    verified_hashes[str(checkpoint.resolve())] = checkpoint_hash
    verified_hashes[str(worker_script.resolve())] = sha256(worker_script)
    return verified_hashes


def artifact_counts(
    campaign: Path, jobs: list[dict[str, Any]], reused_count: int
) -> dict[str, int]:
    """Count completed prediction, normalization, and relaxation artifacts."""
    predictions = normalized = relaxed = 0
    for job in jobs:
        variant, index = job["variant"], int(job["model_index"])
        canonical = (
            campaign / f"predictions/canonical/{variant}/{variant}_model_{index}.pdb"
        )
        confidence = canonical.with_name(f"confidence_{canonical.stem}.json")
        predictions += int(canonical.is_file() and confidence.is_file())
        case = case_path(campaign, job)
        normalized += int(
            (case / "raw/normalized.pdb").is_file()
            and (case / "request.json").is_file()
        )
        process_path = case / "relaxation/process_acceptance.json"
        if (
            process_path.is_file()
            and json.loads(process_path.read_text()).get("status") == "completed"
        ):
            relaxed += 1
    return {
        "predicted_models": predictions,
        "normalized_models": normalized,
        "relaxed_models": relaxed,
        "reused_models": reused_count,
    }
