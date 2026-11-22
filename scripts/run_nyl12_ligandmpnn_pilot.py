#!/usr/bin/env python3
"""Run the frozen, preflight-approved LigandMPNN pilot and audit its outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.nyl12_covalent.ligandmpnn_pilot_audit import (  # noqa: E402
    audit_campaign,
    audit_job,
)

LOGGER = logging.getLogger(__name__)
AUDITOR_SOURCE = Path(sys.modules[audit_job.__module__].__file__).resolve()
REQUIRED_ENVIRONMENT = {
    "CUDA_VISIBLE_DEVICES": "",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
}
EXPECTED_JOBS = 4
EXPECTED_SAMPLES = 24
JOB_TIMEOUT_SECONDS = 1200


def _sha256(path: Path) -> str:
    """Hash an existing file without changing it."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _now() -> str:
    """Return a timezone-aware UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write a deterministic JSON status or report document."""
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _validate_preflight(
    manifest_path: Path, receipt_path: Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate receipt acceptance and all frozen source/input hashes."""
    manifest = json.loads(manifest_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("accepted") is not True:
        raise ValueError("Preflight receipt is not accepted")
    if receipt.get("host") != "violet":
        raise ValueError("Preflight receipt does not authorize host violet")
    if (
        receipt.get("parent_models") != [23]
        or receipt.get("sample_budget") != EXPECTED_SAMPLES
    ):
        raise ValueError("Preflight receipt does not authorize model 23 and 24 samples")
    if _sha256(manifest_path) != receipt.get("manifest_sha256"):
        raise ValueError("Frozen manifest SHA256 does not match preflight receipt")
    source_manifest = Path(receipt["input_source_manifest"])
    if _sha256(source_manifest) != receipt.get("input_source_manifest_sha256"):
        raise ValueError(
            "Input source manifest SHA256 does not match preflight receipt"
        )
    source_data = json.loads(source_manifest.read_text())
    source_root = source_manifest.parent.resolve()
    declared_sources: set[str] = set()
    for relative_path, metadata in source_data["files"].items():
        source_path = (source_root / relative_path).resolve()
        if _sha256(source_path) != metadata["sha256"]:
            raise ValueError(f"Input source SHA256 mismatch: {source_path}")
        declared_sources.add(str(source_path))
    execution_hashes = receipt["execution_sha256"]
    if not execution_hashes:
        raise ValueError("Preflight receipt has no execution source hashes")
    for raw_path, expected in execution_hashes.items():
        path = Path(raw_path)
        if not path.is_absolute():
            raise ValueError(f"Execution hash path must be absolute: {path}")
        if _sha256(path) != expected:
            raise ValueError(f"Execution source SHA256 mismatch: {path}")
    hashed_sources = {str(Path(path).resolve()) for path in execution_hashes}
    required_sources = {
        str(Path(__file__).resolve()),
        str(AUDITOR_SOURCE),
    } | declared_sources
    if not required_sources <= hashed_sources:
        raise ValueError(
            f"Preflight execution hashes omit required sources: "
            f"{sorted(required_sources - hashed_sources)}"
        )
    environment = receipt["environment"]
    if environment != REQUIRED_ENVIRONMENT:
        raise ValueError(
            "Preflight environment overrides do not match the approved values"
        )
    if len(manifest["jobs"]) != EXPECTED_JOBS:
        raise ValueError(f"Frozen manifest must contain {EXPECTED_JOBS} jobs")
    if manifest["config"].get("models") != [23]:
        raise ValueError("Frozen manifest must target parent model 23")
    if manifest.get("expected_sequences_before_deduplication") != EXPECTED_SAMPLES:
        raise ValueError(
            f"Frozen manifest must budget exactly {EXPECTED_SAMPLES} samples"
        )
    for backbone in manifest["backbones"]:
        path = Path(backbone["design_pdb"])
        if _sha256(path) != backbone["design_pdb_sha256"]:
            raise ValueError(f"Prepared backbone SHA256 mismatch: {path}")
    return manifest, receipt


def _job_paths(job: dict[str, Any]) -> tuple[list[str], Path, Path]:
    """Extract the frozen argv, LigandMPNN output folder, and engine cwd."""
    argv = job["argv"]
    try:
        output_index = argv.index("--out_folder")
        output_folder = Path(argv[output_index + 1])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Job {job['job_id']} lacks --out_folder") from exc
    if len(argv) < 2:
        raise ValueError(f"Job {job['job_id']} has an incomplete argv")
    try:
        batch_size = int(argv[argv.index("--batch_size") + 1])
        batches = int(argv[argv.index("--number_of_batches") + 1])
        temperature = float(argv[argv.index("--temperature") + 1])
        seed = int(argv[argv.index("--seed") + 1])
    except (ValueError, IndexError) as exc:
        raise ValueError(
            f"Job {job['job_id']} lacks sampling count/lineage flags"
        ) from exc
    if batch_size * batches != 6:
        raise ValueError(f"Job {job['job_id']} does not request exactly six samples")
    if not 0 < temperature <= 1 or seed < 0:
        raise ValueError(f"Job {job['job_id']} has invalid temperature or seed")
    return argv, output_folder, Path(argv[1]).resolve().parent


def _sample_fasta(samples: list[dict[str, Any]], unique: bool) -> str:
    """Serialize sample sequences with copy strings and deduplication lineage."""
    groups: dict[str, list[dict[str, Any]]] = {}
    for sample in samples:
        sequence = sample["alpha_sequence"] + ":" + sample["beta_sequence"]
        groups.setdefault(sequence, []).append(sample)
    records: list[str] = []
    chosen = groups.values() if unique else [[sample] for sample in samples]
    for group in chosen:
        first = group[0]
        lineage = ",".join(sample["sample_id"] for sample in group)
        header = f">{first['sample_id']} lineage={lineage} parent={first['parent_model']} T={first['temperature']} seed={first['seed']}"
        alpha, beta = first["alpha_sequence"], first["beta_sequence"]
        records.extend([header, ":".join([alpha] * 4 + [beta] * 4)])
    return "\n".join(records) + "\n"


def _write_exports(output: Path, report: dict[str, Any]) -> None:
    """Write campaign JSON, mutation table, and complete/unique FASTA files."""
    _write_json(output / "report.json", report)
    with (output / "samples.csv").open("w", newline="") as handle:
        fields = [
            "sample_id",
            "parent_model",
            "temperature",
            "seed",
            "mutation_count",
            "mutations",
            "alpha_sequence",
            "beta_sequence",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for sample in report["samples"]:
            writer.writerow(
                {
                    key: json.dumps(sample["mutations"], sort_keys=True)
                    if key == "mutations"
                    else sample[key]
                    for key in fields
                }
            )
    (output / "all_samples.fasta").write_text(
        _sample_fasta(report["samples"], unique=False)
    )
    (output / "unique_sequences.fasta").write_text(
        _sample_fasta(report["samples"], unique=True)
    )
    _write_json(
        output / "scientific_acceptance.json",
        {
            "protein_chain_lengths_validated": True,
            "copy_ties_validated": True,
            "protected_identities_validated": True,
            "untied_positions_unchanged": True,
            "catalytic_activity_validated": False,
        },
    )


def run_campaign(
    manifest_path: Path, preflight_path: Path, output: Path
) -> dict[str, Any]:
    """Execute four frozen jobs sequentially, auditing each completed FASTA.

    Args:
        manifest_path: Frozen preparation manifest.
        preflight_path: Accepted receipt binding the manifest, input, code, and environment.
        output: New directory for logs, statuses, and campaign reports.

    Returns:
        Final campaign report containing all 24 audited sequences.

    Raises:
        FileExistsError: The results directory or any frozen job output already exists.
        FileNotFoundError: A preflight-bound input or source file is absent.
        ValueError: Receipt, hashes, job count, sequence budget, or FASTA audit fails.
        subprocess.CalledProcessError: A LigandMPNN job exits unsuccessfully.
        subprocess.TimeoutExpired: A LigandMPNN job exceeds its 20-minute limit.
    """
    manifest, receipt = _validate_preflight(manifest_path, preflight_path)
    if output.exists():
        raise FileExistsError(output)
    job_paths = [_job_paths(job) for job in manifest["jobs"]]
    accepted_pdbs = {
        str(Path(backbone["design_pdb"]).resolve())
        for backbone in manifest["backbones"]
    }
    for job in manifest["jobs"]:
        argv = job["argv"]
        try:
            pdb_path = Path(argv[argv.index("--pdb_path") + 1]).resolve()
        except (ValueError, IndexError) as exc:
            raise ValueError(f"Job {job['job_id']} lacks --pdb_path") from exc
        if str(pdb_path) not in accepted_pdbs:
            raise ValueError(f"Job {job['job_id']} input PDB is not a frozen backbone")
    lineage = [
        (
            float(job["argv"][job["argv"].index("--temperature") + 1]),
            int(job["argv"][job["argv"].index("--seed") + 1]),
        )
        for job in manifest["jobs"]
    ]
    temperatures = {temperature for temperature, _ in lineage}
    seeds = {seed for _, seed in lineage}
    if (
        len(temperatures) != 2
        or len(seeds) != EXPECTED_JOBS
        or len(set(lineage)) != EXPECTED_JOBS
    ):
        raise ValueError(
            "Jobs must be two temperatures crossed with four independent seeds"
        )
    if any(
        sum(temperature == selected for temperature, _ in lineage) != 2
        for selected in temperatures
    ):
        raise ValueError("Each temperature must have exactly two independent seed jobs")
    for job, (_, folder, _) in zip(manifest["jobs"], job_paths, strict=True):
        if folder.exists():
            raise FileExistsError(f"Refusing to overwrite frozen job output: {folder}")
    output.mkdir(parents=True)
    log_dir = output / "jobs"
    log_dir.mkdir()
    status: dict[str, Any] = {
        "status": "running",
        "total_jobs": EXPECTED_JOBS,
        "completed_jobs": 0,
        "audited_samples": 0,
        "manifest_sha256": receipt["manifest_sha256"],
        "started_at": _now(),
        "jobs": [],
    }
    _write_json(output / "status.json", status)
    environment = os.environ.copy()
    environment.update(receipt["environment"])

    try:
        for job, (argv, _, cwd) in zip(manifest["jobs"], job_paths, strict=True):
            job_id = job["job_id"]
            log_path = log_dir / f"{job_id}.log"
            job_status: dict[str, Any] = {
                "job_id": job_id,
                "argv": argv,
                "cwd": str(cwd),
                "log": str(log_path),
                "status": "running",
                "started_at": _now(),
                "returncode": None,
            }
            status["jobs"].append(job_status)
            _write_json(output / "status.json", status)
            try:
                with log_path.open("w") as log:
                    completed = subprocess.run(
                        argv,
                        cwd=cwd,
                        env=environment,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                        timeout=JOB_TIMEOUT_SECONDS,
                    )
                job_status["returncode"] = completed.returncode
                job_status["completed_at"] = _now()
                job_status["fasta_audited_samples"] = len(audit_job(manifest, job))
                job_status["status"] = "completed"
                status["completed_jobs"] += 1
                status["audited_samples"] += job_status["fasta_audited_samples"]
                _write_json(output / "status.json", status)
            except Exception as exc:
                job_status["status"] = "failed"
                job_status["completed_at"] = _now()
                job_status["error"] = f"{type(exc).__name__}: {exc}"
                if isinstance(exc, subprocess.CalledProcessError):
                    job_status["returncode"] = exc.returncode
                status.update(
                    {
                        "status": "failed",
                        "failed_at": _now(),
                        "error": job_status["error"],
                    }
                )
                _write_json(output / "status.json", status)
                _write_json(log_dir / f"{job_id}.status.json", job_status)
                raise
            _write_json(log_dir / f"{job_id}.status.json", job_status)

        report = audit_campaign(manifest_path)
        if report["summary"]["generated_sample_count"] != EXPECTED_SAMPLES:
            raise ValueError(
                "Campaign audit did not return exactly 24 generated samples"
            )
        _write_exports(output, report)
        status.update({"status": "completed", "completed_at": _now()})
        _write_json(output / "status.json", status)
        return report
    except Exception as exc:
        if status["status"] != "failed":
            status.update(
                {
                    "status": "failed",
                    "failed_at": _now(),
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
            _write_json(output / "status.json", status)
        raise


def main() -> None:
    """Parse the frozen manifest/receipt and execute the approved campaign."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run_campaign(
        args.manifest.resolve(), args.preflight.resolve(), args.output.resolve()
    )
    LOGGER.info("Completed frozen 24-sequence LigandMPNN pilot in %s", args.output)


if __name__ == "__main__":
    main()
