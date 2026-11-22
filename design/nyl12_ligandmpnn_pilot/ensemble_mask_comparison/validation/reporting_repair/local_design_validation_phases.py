"""Prediction, normalization, relaxation, and comparison phase handlers."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import queue
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.design_component_normalization import normalize_model  # noqa: E402
from nyl12_covalent.design_geometry_comparison import compare_models  # noqa: E402
from nyl12_covalent.design_geometry_plots import plot_comparison  # noqa: E402
from nyl12_covalent.design_model_relaxation import relax_model  # noqa: E402
from nyl12_covalent.ligandmpnn_pilot_audit import _pdb_sequences  # noqa: E402
from nyl12_covalent.local_design_validation_inputs import (  # noqa: E402
    artifact_counts,
    atomic_json,
    case_path,
    component_path,
    sha256,
    validate_reused,
)

LOGGER = logging.getLogger(__name__)
WORKER_SCRIPT = ROOT / "scripts/run_design_boltz_predictions.py"
BOLTZ_CACHE = Path("/data/jpc/code/boltz")
BOLTZ_PYTHON = Path("/data/jpc/enzyme-substrate-modeling/.pixi/envs/boltz/bin/python")
MAX_SECONDS = 48 * 60 * 60


def _run_worker(
    campaign: Path,
    previous: Path,
    jobs_path: Path,
    status: dict[str, Any],
    status_path: Path,
) -> None:
    """Run the existing local Boltz worker with the frozen validation environment."""
    argv = [
        sys.executable,
        str(WORKER_SCRIPT),
        "--jobs",
        str(jobs_path),
        "--output-root",
        str(campaign / "predictions"),
        "--source-dir",
        str(previous / "source_snapshot"),
        "--cache",
        str(BOLTZ_CACHE),
        "--gpu-devices",
        "0,1",
        "--boltz-python",
        str(BOLTZ_PYTHON),
    ]
    overrides = {
        "PYTHONNOUSERSITE": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "PYTHONPATH": str(previous / "source_snapshot"),
        "CUDA_VISIBLE_DEVICES": "0,1",
    }
    env = os.environ.copy()
    env.update(overrides)
    env_hash = hashlib.sha256(json.dumps(sorted(env.items())).encode()).hexdigest()
    atomic_json(
        campaign / "worker_invocation.json",
        {
            "argv": argv,
            "cwd": str(ROOT),
            "environment_overrides": overrides,
            "effective_environment_sha256": env_hash,
            "stdout_log": str((campaign / "worker.stdout.log").resolve()),
        },
    )
    status.update(
        phase="prediction",
        **artifact_counts(campaign, status["_jobs"], status["reused_models"]),
    )
    atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
    with (campaign / "worker.stdout.log").open("x", encoding="utf-8") as log:
        result = subprocess.run(
            argv,
            cwd=ROOT,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=MAX_SECONDS,
        )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, argv)
    marker = json.loads((campaign / "predictions/WORKER_COMPLETE.json").read_text())
    if marker.get("jobs") != status["expected_new_models"] or marker.get(
        "complete"
    ) != marker.get("jobs"):
        raise RuntimeError(f"Boltz worker completeness mismatch: {marker}")
    status.update(artifact_counts(campaign, status["_jobs"], status["reused_models"]))
    atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})


def _normalize_jobs(
    campaign: Path, jobs: list[dict[str, Any]], config: dict[str, Any]
) -> None:
    """Normalize every new model and publish its content-bound request record."""
    for job in jobs:
        variant, index = job["variant"], int(job["model_index"])
        source = (
            campaign / f"predictions/canonical/{variant}/{variant}_model_{index}.pdb"
        )
        confidence = source.with_name(f"confidence_{source.stem}.json")
        component = component_path(campaign, job)
        if not all(path.is_file() for path in (source, confidence, component)):
            raise FileNotFoundError(
                f"Prediction artifacts incomplete for {variant}/{index}"
            )
        case = case_path(campaign, job)
        case.mkdir(parents=True, exist_ok=False)
        result = normalize_model(source, component, Path(config["mol2"]), case / "raw")
        roles = result["role_mapping"]
        if roles.get("C54") != "C54" or roles.get("O49") != "O49":
            raise ValueError(f"Reactive atom mapping changed for {variant}/{index}")
        normalized = case / "raw/normalized.pdb"
        observed = _pdb_sequences(normalized)
        for chain in "ABCD":
            if observed[chain] != job["expected_alpha_sequence"]:
                raise ValueError(
                    f"Alpha sequence mismatch for {variant}/{index}/{chain}"
                )
        for chain in "EFGH":
            if observed[chain] != job["expected_beta_sequence"]:
                raise ValueError(
                    f"Beta sequence mismatch for {variant}/{index}/{chain}"
                )
        atomic_json(
            case / "request.json",
            {
                "variant": variant,
                "seed": int(job["seed"]),
                "model_index": index,
                "source": str(normalized.resolve()),
                "source_sha256": sha256(normalized),
                "original_prediction": str(source.resolve()),
                "prediction_sha256": sha256(source),
                "component_pickle": str(component.resolve()),
                "component_sha256": sha256(component),
                "normalization_mapping": str((case / "raw/mapping.json").resolve()),
                "normalization_mapping_sha256": sha256(case / "raw/mapping.json"),
                "ligand_count": 4,
            },
        )


def _relax_new_models(
    campaign: Path, jobs: list[dict[str, Any]], config: dict[str, Any]
) -> None:
    """Relax new models with at most two disjoint GPU workers."""
    devices: queue.Queue[str] = queue.Queue()
    for device in ("0", "1"):
        devices.put(device)

    def relax_one(job: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
        device = devices.get()
        try:
            case = case_path(campaign, job)
            acceptance = relax_model(
                case / "raw/normalized.pdb",
                case / "relaxation",
                {**config, "gpu": device},
            )
            process_path = case / "relaxation/process_acceptance.json"
            process = json.loads(process_path.read_text())
            if (
                process.get("status") != "completed"
                or acceptance.get("process_exit_code") != 0
            ):
                raise RuntimeError(
                    f"Relaxation process failed for {job['variant']}/{job['model_index']}"
                )
            if not isinstance(acceptance.get("force_converged"), bool):
                raise ValueError(
                    f"Missing force-convergence result for {job['variant']}/{job['model_index']}"
                )
            return acceptance, process
        finally:
            devices.put(device)

    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(relax_one, job): job for job in jobs}
        for future in as_completed(futures):
            job = futures[future]
            try:
                future.result()
                LOGGER.info(
                    "Relaxation process completed for %s/%s",
                    job["variant"],
                    job["model_index"],
                )
            except Exception as exc:
                failures.append(f"{job['variant']}/{job['model_index']}: {exc!r}")
    if failures:
        raise RuntimeError("Relaxation process failures: " + "; ".join(failures))


def _compare(
    campaign: Path,
    jobs: list[dict[str, Any]],
    reused: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    """Compare raw/relaxed structures and report process/force results separately."""
    entries: list[dict[str, Any]] = []
    acceptance_rows: list[dict[str, Any]] = []
    for job in jobs:
        case = case_path(campaign, job)
        acceptance_path = case / "relaxation/acceptance.json"
        process_path = case / "relaxation/process_acceptance.json"
        acceptance = json.loads(acceptance_path.read_text())
        process = json.loads(process_path.read_text())
        for stage, path in (
            ("raw", case / "raw/normalized.pdb"),
            ("relaxed", case / "relaxation/relaxed.pdb"),
        ):
            entries.append(
                {
                    "variant": job["variant"],
                    "seed": int(job["seed"]),
                    "stage": stage,
                    "path": str(path.resolve()),
                    "ligand_count": 4,
                }
            )
        acceptance_rows.append(
            {
                "variant": job["variant"],
                "model_index": int(job["model_index"]),
                "seed": int(job["seed"]),
                "process_status": process["status"],
                "process_completed": process["status"] == "completed",
                "force_converged": acceptance["force_converged"],
                "final_rms_gradient": acceptance["final_rms_gradient"],
                "acceptance_json": str(acceptance_path.resolve()),
            }
        )
    for item in reused:
        verified = validate_reused(item)
        for stage, path in (("raw", verified["raw"]), ("relaxed", verified["relaxed"])):
            entries.append(
                {
                    "variant": "WT",
                    "seed": int(item["seed"]),
                    "stage": stage,
                    "path": str(path.resolve()),
                    "ligand_count": 4,
                }
            )
        acceptance = verified["acceptance"]
        process = verified["process"]
        acceptance_rows.append(
            {
                "variant": "WT",
                "model_index": int(item["model_index"]),
                "seed": int(item["seed"]),
                "process_status": process["status"],
                "process_completed": process["status"] == "completed",
                "force_converged": acceptance["force_converged"],
                "final_rms_gradient": acceptance["final_rms_gradient"],
                "acceptance_json": str(verified["acceptance_path"].resolve()),
            }
        )
    analysis = campaign / "analysis"
    result = compare_models(entries, Path(config["reference_cif"]), analysis)
    plots = plot_comparison(
        analysis / "occupied_sites.csv",
        analysis / "model_summary.csv",
        analysis / "plots",
        planned_models=len(jobs) + len(reused),
        variant_order=tuple(dict.fromkeys(entry["variant"] for entry in entries)),
    )
    with (analysis / "relaxation_acceptance.csv").open(
        "x", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(acceptance_rows[0]))
        writer.writeheader()
        writer.writerows(acceptance_rows)
    return {
        "geometry_report": result,
        "plots": plots,
        "acceptance_rows": acceptance_rows,
        "entries": entries,
    }
