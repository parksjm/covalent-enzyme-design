#!/usr/bin/env python3
"""Run a local, bounded Boltz and Amber validation campaign."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.local_design_validation_inputs import (  # noqa: E402
    artifact_counts,
    atomic_json,
    sha256,
    validate_jobs,
    validate_reused,
    verify_runtime,
)
from nyl12_covalent.local_design_validation_phases import (  # noqa: E402
    BOLTZ_CACHE,
    BOLTZ_PYTHON,
    WORKER_SCRIPT,
    _compare,
    _normalize_jobs,
    _relax_new_models,
    _run_worker,
)

LOGGER = logging.getLogger(__name__)


def run_validation(campaign: Path, previous_campaign: Path) -> dict[str, Any]:
    """Run a bounded local prediction, normalization, relaxation, and comparison.

    Args:
        campaign: New validation directory with prepared inputs and configuration.
        previous_campaign: Completed validation campaign providing WT reuse and
            frozen Boltz runtime provenance.

    Returns:
        Completion metadata with source hashes and model-level counts.

    Raises:
        FileExistsError: If status/completion or result directories already exist.
        ValueError: If manifest lineage, reuse hashes, or inputs are invalid.
        subprocess.CalledProcessError: If the prediction worker fails.
        RuntimeError: If normalization, relaxation, or comparison fails.
    """
    campaign, previous = campaign.resolve(), previous_campaign.resolve()
    status_path, completion_path = (
        campaign / "status.json",
        campaign / "completion.json",
    )
    if status_path.exists() or completion_path.exists():
        raise FileExistsError("Existing status/completion requires explicit review")
    for name in (
        "predictions",
        "models",
        "analysis",
        "worker_invocation.json",
        "worker.stdout.log",
    ):
        if (campaign / name).exists():
            raise FileExistsError(campaign / name)
    status: dict[str, Any] = {
        "status": "running",
        "phase": "preflight",
        "predicted_models": 0,
        "normalized_models": 0,
        "relaxed_models": 0,
        "reused_models": 0,
        "expected_new_models": 0,
        "expected_reused_models": 3,
        "expected_total_models": 0,
    }
    atomic_json(status_path, status)
    try:
        manifest_path = campaign / "prepared/candidate_manifest.json"
        jobs_path = campaign / "prepared/jobs.csv"
        config_path = campaign / "configuration/relaxation.json"
        manifest, config = (
            json.loads(manifest_path.read_text()),
            json.loads(config_path.read_text()),
        )
        jobs = validate_jobs(campaign, manifest, ROOT)
        reused = manifest["reused_models"]
        verified_reuse = [validate_reused(item) for item in reused]
        runtime_hashes = verify_runtime(
            previous, WORKER_SCRIPT, BOLTZ_CACHE, BOLTZ_PYTHON
        )
        for key in ("release", "mol2", "mapping", "reference_cif"):
            source = Path(config[key])
            valid_type = source.is_dir() if key == "release" else source.is_file()
            if not source.is_absolute() or not valid_type:
                raise FileNotFoundError(source)
        input_hashes = {
            str(path.resolve()): sha256(path)
            for path in (
                manifest_path,
                jobs_path,
                config_path,
                previous / "status.json",
                previous / "completion.json",
            )
        }
        source_hashes = {**runtime_hashes, **input_hashes}
        for item, verified in zip(reused, verified_reuse, strict=True):
            source_hashes[str(verified["raw"].resolve())] = verified["raw_sha256"]
            source_hashes[str(verified["relaxed"].resolve())] = sha256(
                verified["relaxed"]
            )
            source_hashes[str(verified["acceptance_path"].resolve())] = sha256(
                verified["acceptance_path"]
            )
        source_hashes[str(Path(__file__).resolve())] = sha256(Path(__file__).resolve())
        status.update(
            expected_new_models=len(jobs),
            expected_total_models=len(jobs) + len(reused),
            reused_models=len(reused),
            _jobs=jobs,
        )
        atomic_json(campaign / "source_hashes.json", source_hashes)
        atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        _run_worker(campaign, previous, jobs_path, status, status_path)
        status["phase"] = "normalization"
        atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        _normalize_jobs(campaign, jobs, config)
        status.update(
            phase="normalization", **artifact_counts(campaign, jobs, len(reused))
        )
        atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        if status["normalized_models"] != len(jobs):
            raise RuntimeError("Not all new models normalized")
        status["phase"] = "relaxation"
        atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        _relax_new_models(campaign, jobs, config)
        status.update(**artifact_counts(campaign, jobs, len(reused)))
        atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        if status["relaxed_models"] != len(jobs):
            raise RuntimeError("Not all new relaxation processes completed")
        status["phase"] = "comparison"
        atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        comparison = _compare(campaign, jobs, reused, config)
        completion = {
            "status": "completed",
            "new_models": len(jobs),
            "reused_models": len(reused),
            "total_model_seeds": len(jobs) + len(reused),
            "predicted_models": status["predicted_models"],
            "normalized_models": status["normalized_models"],
            "relaxed_models": status["relaxed_models"],
            "process_completed_models": sum(
                row["process_completed"] for row in comparison["acceptance_rows"]
            ),
            "force_converged_models": sum(
                row["force_converged"] for row in comparison["acceptance_rows"]
            ),
            "sample_unit": "Model seed; four catalytic sites are nested within each model.",
            "interpretation": "Geometry and minimization do not validate catalytic activity.",
            "source_hashes": str((campaign / "source_hashes.json").resolve()),
            "geometry_report": comparison["geometry_report"],
            "plots": comparison["plots"],
        }
        atomic_json(completion_path, completion)
        status.update(
            status="completed",
            phase="complete",
            **artifact_counts(campaign, jobs, len(reused)),
        )
        atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        return completion
    except Exception as exc:
        status.update(status="failed", error=repr(exc), error_type=type(exc).__name__)
        try:
            if "_jobs" in status:
                status.update(
                    artifact_counts(campaign, status["_jobs"], status["reused_models"])
                )
        finally:
            atomic_json(status_path, {k: v for k, v in status.items() if k != "_jobs"})
        raise


def main() -> None:
    """Run the local design validation controller from prepared campaign files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--previous-campaign", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    run_validation(args.campaign, args.previous_campaign)


if __name__ == "__main__":
    main()
