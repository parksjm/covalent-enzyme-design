"""Tests for the bounded local design validation controller."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import threading
import time
from unittest.mock import patch

import pytest

from scripts import run_local_design_validation as control
from nyl12_covalent import local_design_validation_phases as phases
from nyl12_covalent import local_design_validation_inputs as validation_inputs

ALPHA = "A" * 232
BETA = "T" * 96


def _reused_models(base: Path) -> list[dict]:
    """Create three immutable WT reuse records with separate process evidence."""
    records = []
    campaign_root = base / "boltz_validation"
    for index in (0, 2, 4):
        seed = index + 1
        folder = campaign_root / "models/WT" / f"seed_{seed:02d}"
        raw_folder = folder / "raw"
        raw_folder.mkdir(parents=True)
        raw = raw_folder / "normalized.pdb"
        raw.write_text(f"raw-{index}")
        mapping = raw_folder / "mapping.json"
        mapping.write_text(f"mapping-{index}")
        yaml = campaign_root / "prepared/inputs/WT.yaml"
        yaml.parent.mkdir(parents=True, exist_ok=True)
        yaml.write_text("WT YAML")
        source_copy = folder / "relaxation/inputs/source_model.pdb"
        source_copy.parent.mkdir(parents=True)
        source_copy.write_bytes(raw.read_bytes())
        relaxed = folder / "relaxation/relaxed.pdb"
        relaxed.parent.mkdir(parents=True, exist_ok=True)
        relaxed.write_text(f"relaxed-{index}")
        acceptance = folder / "relaxation/acceptance.json"
        acceptance.write_text(
            json.dumps(
                {
                    "endpoint": str(relaxed),
                    "source": str(source_copy),
                    "process_exit_code": 0,
                    "force_converged": True,
                    "final_rms_gradient": 0.01,
                }
            )
        )
        (folder / "relaxation/process_acceptance.json").write_text(
            json.dumps({"status": "completed"})
        )
        records.append(
            {
                "variant": "WT",
                "model_index": index,
                "seed": 900 + index,
                "raw_path": str(raw),
                "relaxed_path": str(relaxed),
                "acceptance_json": str(acceptance),
                "source_sha256": {
                    "raw": hashlib.sha256(raw.read_bytes()).hexdigest(),
                    "relaxed": hashlib.sha256(relaxed.read_bytes()).hexdigest(),
                    "acceptance": hashlib.sha256(acceptance.read_bytes()).hexdigest(),
                    "mapping": hashlib.sha256(mapping.read_bytes()).hexdigest(),
                    "yaml": hashlib.sha256(yaml.read_bytes()).hexdigest(),
                },
            }
        )
    return records


def _campaign(tmp_path: Path, n_jobs: int = 3) -> tuple[Path, Path, list[dict]]:
    """Prepare a minimal campaign contract with two new WT and one design job."""
    campaign = tmp_path / "validation"
    previous = tmp_path / "previous"
    for folder in (
        campaign / "prepared",
        campaign / "configuration",
        previous / "source_snapshot",
    ):
        folder.mkdir(parents=True, exist_ok=True)
    records = _reused_models(previous)
    indices = [("WT", 1), ("WT", 3), ("variant_01", 0)][:n_jobs]
    jobs = []
    for serial, (variant, index) in enumerate(indices):
        yaml = campaign / "prepared" / f"{variant}_{index}.yaml"
        yaml.write_text("input")
        jobs.append(
            {
                "variant": variant,
                "model_index": index,
                "seed": 1000 + serial,
                "yaml_path": str(yaml),
                "expected_alpha_sequence": ALPHA,
                "expected_beta_sequence": BETA,
            }
        )
    (campaign / "prepared/candidate_manifest.json").write_text(
        json.dumps({"jobs": jobs, "reused_models": records})
    )
    with (campaign / "prepared/jobs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "campaign",
                "sequence_name",
                "model_index",
                "seed",
                "yaml_path",
            ],
        )
        writer.writeheader()
        for job in jobs:
            writer.writerow(
                {
                    "campaign": "validation",
                    "sequence_name": job["variant"],
                    **{key: job[key] for key in ("model_index", "seed", "yaml_path")},
                }
            )
    for name in ("release", "mol2", "mapping", "reference_cif"):
        path = campaign / "configuration" / name
        if name == "release":
            path.mkdir()
        else:
            path.write_text(name)
    config = {
        "release": str(campaign / "configuration/release"),
        "mol2": str(campaign / "configuration/mol2"),
        "mapping": str(campaign / "configuration/mapping"),
        "reference_cif": str(campaign / "configuration/reference_cif"),
    }
    (campaign / "configuration/relaxation.json").write_text(json.dumps(config))
    (previous / "completion.json").write_text(json.dumps({"process_complete": True}))
    (previous / "status.json").write_text(json.dumps({"status": "completed"}))
    (previous / "source_snapshot/run_boltz_cli.py").write_text("frozen runtime")
    return campaign, previous, jobs


def _fake_worker(
    campaign: Path, previous: Path, jobs_path: Path, status: dict, status_path: Path
) -> None:
    """Create all canonical prediction and component artifacts without Boltz."""
    for job in status["_jobs"]:
        variant, index = job["variant"], job["model_index"]
        canonical = (
            campaign / f"predictions/canonical/{variant}/{variant}_model_{index}.pdb"
        )
        canonical.parent.mkdir(parents=True, exist_ok=True)
        canonical.write_text(f"pred-{variant}-{index}")
        canonical.with_name(f"confidence_{canonical.stem}.json").write_text("{}")
        component = validation_inputs.component_path(campaign, job)
        component.parent.mkdir(parents=True, exist_ok=True)
        component.write_text("trusted component")
    prediction = campaign / "predictions/WORKER_COMPLETE.json"
    prediction.parent.mkdir(parents=True, exist_ok=True)
    prediction.write_text(
        json.dumps({"jobs": len(status["_jobs"]), "complete": len(status["_jobs"])})
    )
    status.update(
        control.artifact_counts(campaign, status["_jobs"], status["reused_models"])
    )
    control.atomic_json(
        status_path, {key: value for key, value in status.items() if key != "_jobs"}
    )


def _fake_normalize(source: Path, component: Path, mol2: Path, output: Path) -> dict:
    """Preserve source bytes while emitting a graph-map stub."""
    output.mkdir(parents=True)
    (output / "normalized.pdb").write_bytes(source.read_bytes())
    (output / "mapping.json").write_text("{}")
    return {"role_mapping": {"C54": "C54", "O49": "O49"}}


def _fake_relax(source: Path, output: Path, config: dict) -> dict:
    """Emit process-complete acceptance, including one unconverged example."""
    output.mkdir(parents=True)
    endpoint = output / "relaxed.pdb"
    endpoint.write_text(source.read_text())
    converged = "variant_01" not in str(source)
    acceptance = {
        "endpoint": str(endpoint),
        "source": str(source),
        "process_exit_code": 0,
        "force_converged": converged,
        "final_rms_gradient": 0.01 if converged else 0.08,
    }
    (output / "acceptance.json").write_text(json.dumps(acceptance))
    (output / "process_acceptance.json").write_text(json.dumps({"status": "completed"}))
    return acceptance


def _fake_compare(entries: list[dict], reference: Path, output: Path) -> dict:
    """Check both stages and all expected model lineages without geometry tools."""
    assert len(entries) == 12
    assert sum(entry["stage"] == "raw" for entry in entries) == 6
    output.mkdir()
    return {"model_count": len(entries) // 2}


def _run_success(campaign: Path, previous: Path, wrong_sequence: bool = False) -> dict:
    """Run all controller phases with engine entry points mocked."""
    observed = {**dict.fromkeys("ABCD", ALPHA), **dict.fromkeys("EFGH", BETA)}
    if wrong_sequence:
        observed["A"] = "C" * 232
    with (
        patch.object(control, "_run_worker", side_effect=_fake_worker),
        patch.object(control, "verify_runtime", return_value={"runtime": "hash"}),
        patch.object(phases, "normalize_model", side_effect=_fake_normalize),
        patch.object(phases, "_pdb_sequences", return_value=observed),
        patch.object(phases, "relax_model", side_effect=_fake_relax),
        patch.object(phases, "compare_models", side_effect=_fake_compare),
        patch.object(phases, "plot_comparison", return_value={"plotted": True}),
    ):
        return control.run_validation(campaign, previous)


def test_local_campaign_runs_and_keeps_force_result_separate(tmp_path: Path) -> None:
    """All model seeds compare even when a completed minimization is unconverged."""
    campaign, previous, _ = _campaign(tmp_path)

    completion = _run_success(campaign, previous)

    assert completion["new_models"] == 3
    assert completion["reused_models"] == 3
    assert completion["process_completed_models"] == 6
    assert completion["force_converged_models"] == 5
    assert completion["sample_unit"].startswith("Model seed")
    rows = list(
        csv.DictReader((campaign / "analysis/relaxation_acceptance.csv").open())
    )
    assert len(rows) == 6
    variant_row = next(row for row in rows if row["variant"] == "variant_01")
    assert variant_row["process_completed"] == "True"
    assert variant_row["force_converged"] == "False"
    assert json.loads((campaign / "status.json").read_text())["status"] == "completed"


def test_worker_command_and_environment_are_frozen(tmp_path: Path) -> None:
    """Prediction subprocess uses the approved script, source snapshot, and GPUs."""
    campaign, previous, jobs = _campaign(tmp_path)
    status = {"_jobs": jobs, "reused_models": 3, "expected_new_models": len(jobs)}
    status_path = campaign / "status.json"
    status_path.write_text("{}")
    jobs_path = campaign / "prepared/jobs.csv"

    def complete(argv: list[str], **kwargs: object):
        marker = campaign / "predictions/WORKER_COMPLETE.json"
        marker.parent.mkdir(parents=True)
        marker.write_text(json.dumps({"jobs": len(jobs), "complete": len(jobs)}))
        return type("Result", (), {"returncode": 0})()

    with patch.object(phases.subprocess, "run", side_effect=complete) as run:
        control._run_worker(campaign, previous, jobs_path, status, status_path)

    argv = run.call_args.args[0]
    env = run.call_args.kwargs["env"]
    assert argv[0] == control.sys.executable
    assert argv[1] == str(control.WORKER_SCRIPT)
    assert "0,1" in argv
    assert env["PYTHONPATH"] == str(previous / "source_snapshot")
    assert env["CUDA_VISIBLE_DEVICES"] == "0,1"
    assert run.call_args.kwargs["timeout"] == phases.MAX_SECONDS
    record = json.loads((campaign / "worker_invocation.json").read_text())
    assert record["effective_environment_sha256"]


def test_wrong_sequence_is_rejected_and_status_records_failure(tmp_path: Path) -> None:
    """Identity errors fail the phase and leave an atomic failed status."""
    campaign, previous, _ = _campaign(tmp_path)

    with pytest.raises(ValueError, match="Alpha sequence"):
        _run_success(campaign, previous, wrong_sequence=True)

    status = json.loads((campaign / "status.json").read_text())
    assert status["status"] == "failed"
    assert status["phase"] == "normalization"
    assert "error" in status


def test_reused_source_hash_mismatch_blocks_campaign(tmp_path: Path) -> None:
    """A changed reused WT raw source is rejected before any new prediction."""
    campaign, previous, _ = _campaign(tmp_path)
    manifest_path = campaign / "prepared/candidate_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["reused_models"][0]["source_sha256"]["raw"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="SHA256 mismatch"):
        _run_success(campaign, previous)

    status = json.loads((campaign / "status.json").read_text())
    assert status["status"] == "failed"
    assert status["phase"] == "preflight"


def test_worker_failure_is_recorded_without_retry(tmp_path: Path) -> None:
    """A failed prediction worker leaves failure state and no silent restart."""
    campaign, previous, _ = _campaign(tmp_path)
    with (
        patch.object(control, "verify_runtime", return_value={"runtime": "hash"}),
        patch.object(control, "_run_worker", side_effect=RuntimeError("worker failed")),
        pytest.raises(RuntimeError, match="worker failed"),
    ):
        control.run_validation(campaign, previous)
    assert json.loads((campaign / "status.json").read_text())["status"] == "failed"
    with pytest.raises(FileExistsError, match="status/completion"):
        control.run_validation(campaign, previous)


def test_relaxation_workers_own_two_gpus_and_finish_all_jobs(tmp_path: Path) -> None:
    """The relaxation queue never assigns a GPU to two simultaneous jobs."""
    campaign, _, jobs = _campaign(tmp_path)
    for job in jobs:
        case = validation_inputs.case_path(campaign, job)
        (case / "raw").mkdir(parents=True)
        (case / "raw/normalized.pdb").write_text("source")
    lock = threading.Lock()
    active: set[str] = set()
    seen: set[str] = set()
    maximum = 0

    def tracked(source: Path, output: Path, config: dict):
        nonlocal maximum
        gpu = config["gpu"]
        with lock:
            assert gpu not in active
            active.add(gpu)
            seen.add(gpu)
            maximum = max(maximum, len(active))
        time.sleep(0.02)
        result = _fake_relax(source, output, config)
        with lock:
            active.remove(gpu)
        return result

    with patch.object(phases, "relax_model", side_effect=tracked):
        control._relax_new_models(campaign, jobs, {})

    assert seen == {"0", "1"}
    assert maximum <= 2


def test_jobs_csv_must_match_manifest_lineage(tmp_path: Path) -> None:
    """A changed seed or YAML mapping cannot reach the worker manifest."""
    campaign, _, jobs = _campaign(tmp_path)
    manifest = json.loads((campaign / "prepared/candidate_manifest.json").read_text())
    jobs_path = campaign / "prepared/jobs.csv"
    rows = list(csv.DictReader(jobs_path.open()))
    rows[0]["seed"] = str(int(rows[0]["seed"]) + 1)
    with jobs_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

    with pytest.raises(ValueError, match="lineage mismatch"):
        validation_inputs.validate_jobs(campaign, manifest, control.ROOT)
    assert len(jobs) == 3


def test_runtime_snapshot_hashes_are_checked(tmp_path: Path) -> None:
    """Native runtime, coordinator, and checkpoint receipts bind real inputs."""
    previous = tmp_path / "previous"
    snapshot = previous / "source_snapshot"
    checked_source = snapshot / "checked_runtime/module.py"
    generation_source = snapshot / "scripts/run_design_boltz_predictions.py"
    snapshot.mkdir(parents=True)
    checked_source.parent.mkdir()
    generation_source.parent.mkdir()
    (snapshot / "run_boltz_cli.py").write_text("boltz runtime")
    checked_source.write_text("frozen checked source")
    generation_source.write_text("generation source")
    worker = tmp_path / "validation_snapshot/scripts/run_design_boltz_predictions.py"
    worker.parent.mkdir(parents=True)
    worker.write_text("worker")
    cache = tmp_path / "cache"
    cache.mkdir()
    python = tmp_path / "python"
    python.write_text("python")
    checkpoint = cache / "boltz2_conf.ckpt"
    checkpoint.write_text("qualified checkpoint")
    candidate = previous / "prepared/candidate_manifest.json"
    candidate.parent.mkdir()
    candidate.write_text("{}")
    config = previous / "configuration/relaxation.json"
    config.parent.mkdir()
    config.write_text("{}")
    checked_manifest = previous / "dispatch/checked_runtime_manifest.json"
    checked_manifest.parent.mkdir()
    checked_manifest.write_text(
        json.dumps(
            {
                "sha256": {
                    "source_snapshot/checked_runtime/module.py": hashlib.sha256(
                        checked_source.read_bytes()
                    ).hexdigest()
                }
            }
        )
    )
    plan = previous / "dispatch/plan.json"
    plan.write_text(
        json.dumps(
            {
                "source_sha256": {
                    str(generation_source): hashlib.sha256(
                        generation_source.read_bytes()
                    ).hexdigest()
                }
            }
        )
    )
    runtime = previous / "dispatch/runtime_versions.json"
    runtime.write_text(
        json.dumps(
            {
                "same_boltz_source": True,
                "boltz_distribution": "2.2.1",
                "checkpoint_sha256": hashlib.sha256(
                    checkpoint.read_bytes()
                ).hexdigest(),
                "checkpoint_bytes": checkpoint.stat().st_size,
                "violet": {"torch": "2.9.1+cu128"},
            }
        )
    )
    coordinator = previous / "jobs/coordinator_checked/status.json"
    coordinator.parent.mkdir(parents=True)
    input_hashes = {
        str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (checked_manifest, candidate, config)
    }
    coordinator.write_text(
        json.dumps(
            {
                "exit_status": 0,
                "command": ["coordinator", "--campaign", str(previous.resolve())],
                "input_sha256": input_hashes,
            }
        )
    )
    campaign_status = previous / "status.json"
    campaign_status.write_text(
        json.dumps(
            {
                "status": "completed",
                "planned_models": 45,
                "normalized_models": 45,
                "relaxed_models": 45,
            }
        )
    )
    receipt = previous / "verification/checkpoint_receipt.json"
    receipt.parent.mkdir()
    receipt.write_text(
        json.dumps(
            {
                "prediction_runtime": {
                    "source_sha256": hashlib.sha256(runtime.read_bytes()).hexdigest(),
                    "checkpoint_sha256": hashlib.sha256(
                        checkpoint.read_bytes()
                    ).hexdigest(),
                }
            }
        )
    )

    checked = validation_inputs.verify_runtime(previous, worker, cache, python)

    assert (
        checked[str(checked_source.resolve())]
        == hashlib.sha256(checked_source.read_bytes()).hexdigest()
    )
    checked_source.write_text("modified source")
    with pytest.raises(ValueError, match="hash mismatch"):
        validation_inputs.verify_runtime(previous, worker, cache, python)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ("not_wt", "Only WT models"),
        ("relaxed_hash", "relaxed source SHA256 mismatch"),
        ("acceptance_hash", "acceptance source SHA256 mismatch"),
        ("mapping_hash", "mapping source SHA256 mismatch"),
        ("yaml_hash", "yaml source SHA256 mismatch"),
        ("endpoint", "endpoint path mismatch"),
        ("source_copy", "does not match raw model"),
        ("process", "did not complete"),
        ("force", "force-convergence result"),
    ],
)
def test_reused_acceptance_evidence_is_required(
    tmp_path: Path, change: str, message: str
) -> None:
    """Reuse cannot proceed with altered source, endpoint, process, or force data."""
    record = _reused_models(tmp_path)[0]
    acceptance_path = Path(record["acceptance_json"])
    acceptance = json.loads(acceptance_path.read_text())
    if change == "not_wt":
        record["variant"] = "variant_01"
    elif change.endswith("_hash"):
        source_key = change.removesuffix("_hash")
        record["source_sha256"][source_key] = "0" * 64
    elif change == "endpoint":
        acceptance["endpoint"] = str(tmp_path / "other.pdb")
        acceptance_path.write_text(json.dumps(acceptance))
    elif change == "source_copy":
        Path(acceptance["source"]).write_text("altered copy")
    elif change == "process":
        acceptance_path.with_name("process_acceptance.json").write_text(
            json.dumps({"status": "failed"})
        )
    elif change == "force":
        acceptance.pop("force_converged")
        acceptance_path.write_text(json.dumps(acceptance))
    if change in {"endpoint", "force"}:
        record["source_sha256"]["acceptance"] = hashlib.sha256(
            acceptance_path.read_bytes()
        ).hexdigest()

    with pytest.raises(ValueError, match=message):
        validation_inputs.validate_reused(record)


def test_reused_hash_record_must_be_a_complete_mapping(tmp_path: Path) -> None:
    """A scalar digest cannot silently replace the five named source hashes."""
    record = _reused_models(tmp_path)[0]
    record["source_sha256"] = "0" * 64

    with pytest.raises(
        ValueError, match="must map raw/relaxed/acceptance/mapping/yaml"
    ):
        validation_inputs.validate_reused(record)
