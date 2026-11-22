"""Execute and compare frozen LigandMPNN design-mask jobs."""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import logging
import math
import os
from datetime import datetime, timezone
from pathlib import Path
import re
import subprocess
from typing import Any

from .ligandmpnn_pilot_audit import _pdb_sequences, audit_job

LOGGER = logging.getLogger(__name__)
EXPECTED_JOBS = 8
EXPECTED_SAMPLES = 48
SAMPLES_PER_JOB = 6
JOB_TIMEOUT_SECONDS = 1200
REQUIRED_ENVIRONMENT = {
    "CUDA_VISIBLE_DEVICES": "",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, data: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _flag(argv: list[str], flag: str, job_id: str) -> str:
    try:
        return argv[argv.index(flag) + 1]
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Job {job_id} is missing {flag}") from exc


def _validate_manifest(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    if manifest.get("schema_version") != 1:
        raise ValueError("Manifest schema_version must be 1")
    if manifest.get("expected_jobs") != EXPECTED_JOBS:
        raise ValueError(f"Manifest expected_jobs must be {EXPECTED_JOBS}")
    if manifest.get("expected_samples") != EXPECTED_SAMPLES:
        raise ValueError(f"Manifest expected_samples must be {EXPECTED_SAMPLES}")
    if manifest.get("environment") != REQUIRED_ENVIRONMENT:
        raise ValueError("Manifest environment does not match the frozen CPU settings")
    jobs = manifest.get("jobs")
    if not isinstance(jobs, list) or len(jobs) != EXPECTED_JOBS:
        raise ValueError(f"Manifest must contain exactly {EXPECTED_JOBS} jobs")
    hashes = manifest.get("source_sha256")
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError("Manifest source_sha256 must be a nonempty path/hash map")
    for raw_path, expected in hashes.items():
        path = Path(raw_path)
        if not path.is_absolute():
            raise ValueError(f"Source hash path must be absolute: {path}")
        if not path.is_file():
            raise FileNotFoundError(path)
        if _sha256(path) != expected:
            raise ValueError(f"Source SHA256 mismatch: {path}")

    job_ids: set[str] = set()
    output_paths: set[Path] = set()
    combinations: set[tuple[str, str, float]] = set()
    parent_values: set[str] = set()
    arm_values: set[str] = set()
    temperatures: set[float] = set()
    seeds_by_parent_temperature: dict[tuple[str, float], dict[str, int]] = {}
    sequence_families: tuple[str, str] | None = None
    for job in jobs:
        job_id = job.get("job_id")
        parent, arm = job.get("parent"), job.get("arm")
        if not isinstance(job_id, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", job_id):
            raise ValueError("Each job requires a filesystem-safe job_id")
        if job_id in job_ids:
            raise ValueError(f"Duplicate job_id: {job_id}")
        job_ids.add(job_id)
        if not isinstance(parent, (str, int)) or not str(parent):
            raise ValueError(f"Job {job_id} requires a parent identifier")
        if not isinstance(arm, str) or not arm:
            raise ValueError(f"Job {job_id} requires a nonempty arm")
        parent, arm = str(parent), str(arm)
        parent_values.add(parent)
        arm_values.add(arm)
        argv = job.get("argv")
        if (
            not isinstance(argv, list)
            or len(argv) < 2
            or not all(isinstance(value, str) for value in argv)
        ):
            raise ValueError(f"Job {job_id} argv must be a string argument list")
        out_folder = Path(_flag(argv, "--out_folder", job_id)).resolve()
        if any(
            out_folder == existing
            or out_folder in existing.parents
            or existing in out_folder.parents
            for existing in output_paths
        ):
            raise ValueError(f"Overlapping job output folder: {out_folder}")
        output_paths.add(out_folder)
        if out_folder.exists():
            raise FileExistsError(out_folder)
        if (
            int(_flag(argv, "--batch_size", job_id))
            * int(_flag(argv, "--number_of_batches", job_id))
            != SAMPLES_PER_JOB
        ):
            raise ValueError(f"Job {job_id} must request exactly six samples")
        temperature = float(_flag(argv, "--temperature", job_id))
        if not math.isfinite(temperature) or temperature <= 0:
            raise ValueError(f"Job {job_id} temperature must be positive")
        seed = int(_flag(argv, "--seed", job_id))
        if seed <= 0:
            raise ValueError(f"Job {job_id} seed must be positive")
        combination = (parent, arm, temperature)
        if combination in combinations:
            raise ValueError(f"Duplicate parent/arm/temperature: {combination}")
        combinations.add(combination)
        temperatures.add(temperature)
        seeds_by_parent_temperature.setdefault((parent, temperature), {})[arm] = seed
        input_pdb = Path(_flag(argv, "--pdb_path", job_id)).resolve()
        checkpoint = Path(_flag(argv, "--checkpoint_ligand_mpnn", job_id)).resolve()
        engine = Path(argv[1]).resolve()
        required = {str(input_pdb), str(checkpoint), str(engine)}
        if not required <= hashes.keys():
            raise ValueError(
                f"Job {job_id} source hashes omit PDB/checkpoint/engine: "
                f"{sorted(required - hashes.keys())}"
            )
        sequences = _pdb_sequences(input_pdb)
        if any(sequences[chain] != sequences["A"] for chain in "BCD"):
            raise ValueError(f"Job {job_id} alpha chains A-D differ")
        if any(sequences[chain] != sequences["E"] for chain in "FGH"):
            raise ValueError(f"Job {job_id} beta chains E-H differ")
        families = sequences["A"], sequences["E"]
        if sequence_families is None:
            sequence_families = families
        elif sequence_families != families:
            raise ValueError(
                f"Job {job_id} PDB family sequences differ from other jobs"
            )
        audit_manifest = job.get("audit_manifest")
        if not isinstance(audit_manifest, dict):
            raise ValueError(f"Job {job_id} requires audit_manifest")
        mask, config = audit_manifest.get("design_mask"), audit_manifest.get("config")
        if not isinstance(mask, dict) or not isinstance(config, dict):
            raise ValueError(f"Job {job_id} audit_manifest lacks mask/config")
        for key in ("alpha_positions", "beta_positions"):
            if not isinstance(mask.get(key), list):
                raise ValueError(f"Job {job_id} design_mask lacks {key}")
        for key in ("protected_alpha_positions", "protected_beta_positions"):
            if not isinstance(config.get(key), list):
                raise ValueError(f"Job {job_id} config lacks {key}")
    if len(parent_values) != 2 or len(arm_values) != 2 or len(temperatures) != 2:
        raise ValueError("Eight jobs must cover two parents, arms, and temperatures")
    if len(combinations) != len(parent_values) * len(arm_values) * len(temperatures):
        raise ValueError("Manifest must contain every parent/arm/temperature pair")
    if any(
        set(arms) != arm_values or len(set(arms.values())) != 1
        for arms in seeds_by_parent_temperature.values()
    ):
        raise ValueError("Each parent/temperature seed must match across both arms")
    if (
        len(
            {next(iter(arms.values())) for arms in seeds_by_parent_temperature.values()}
        )
        != 4
    ):
        raise ValueError("Each parent/temperature stratum must use its own seed")
    return jobs


def _run_job(job: dict[str, Any], output: Path) -> dict[str, Any]:
    job_id, argv = job["job_id"], job["argv"]
    process_dir = output / "jobs" / job_id
    process_dir.mkdir(parents=True)
    environment = os.environ.copy()
    environment.update(REQUIRED_ENVIRONMENT)
    started = _now()
    status: dict[str, Any] = {
        "job_id": job_id,
        "parent": job["parent"],
        "arm": job["arm"],
        "argv": argv,
        "cwd": str(Path(argv[1]).resolve().parent),
        "started_utc": started,
    }
    try:
        result = subprocess.run(
            argv,
            cwd=status["cwd"],
            env=environment,
            capture_output=True,
            text=True,
            timeout=JOB_TIMEOUT_SECONDS,
            check=False,
        )
        (process_dir / "stdout.log").write_text(result.stdout)
        (process_dir / "stderr.log").write_text(result.stderr)
        status["exit_code"] = result.returncode
        if result.returncode:
            status["status"] = "failed"
    except (OSError, subprocess.TimeoutExpired) as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        (process_dir / "stdout.log").write_text(
            stdout.decode() if isinstance(stdout, bytes) else stdout
        )
        (process_dir / "stderr.log").write_text(
            stderr.decode() if isinstance(stderr, bytes) else stderr
        )
        status.update(status="failed", error=repr(exc))
        status["ended_utc"] = _now()
        _write_json(process_dir / "process.json", status)
        raise
    status["ended_utc"] = _now()
    _write_json(process_dir / "process.json", status)
    if status.get("status") == "failed":
        raise subprocess.CalledProcessError(result.returncode, argv)
    return status


def _export(output: Path, report: dict[str, Any]) -> None:
    _write_json(output / "report.json", report)
    fields = [
        "sample_id",
        "job_id",
        "parent_model",
        "parent",
        "arm",
        "temperature",
        "seed",
        "header",
        "mutation_count",
        "mutations",
        "alpha_sequence",
        "beta_sequence",
    ]
    with (output / "samples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for sample in report["samples"]:
            writer.writerow(
                {
                    **sample,
                    "mutations": json.dumps(sample["mutations"], sort_keys=True),
                }
            )
    groups: dict[str, list[dict[str, Any]]] = {}
    for sample in report["samples"]:
        key = sample["alpha_sequence"] + ":" + sample["beta_sequence"]
        groups.setdefault(key, []).append(sample)
    for filename, entries in (
        ("all_samples.fasta", [[sample] for sample in report["samples"]]),
        ("unique_sequences.fasta", list(groups.values())),
    ):
        lines: list[str] = []
        for lineage in entries:
            first = lineage[0]
            ids = ",".join(item["sample_id"] for item in lineage)
            header = (
                f">{first['sample_id']} lineage={ids} "
                f"parent={first['parent']} arm={first['arm']}"
            )
            sequence = ":".join(
                [first["alpha_sequence"]] * 4 + [first["beta_sequence"]] * 4
            )
            lines.extend((header, sequence))
        (output / filename).write_text("\n".join(lines) + "\n")


def run_comparison(
    manifest_path: Path, output: Path, workers: int = 2
) -> dict[str, Any]:
    """Run and audit the frozen eight-job, two-arm design-mask comparison.

    Args:
        manifest_path: Version-1 manifest with source hashes and eight jobs.
        output: New results directory; existing paths are refused.
        workers: CPU job concurrency, limited to one or two workers.

    Returns:
        Audited samples, per-arm counts, job statuses, and sequence uniqueness.

    Raises:
        FileExistsError: Output or any job output directory already exists.
        FileNotFoundError: A declared source file is absent.
        ValueError: Manifest, source hashes, PDB families, or audited outputs fail.
        subprocess.CalledProcessError: A LigandMPNN job exits unsuccessfully.
    """
    if workers not in (1, 2):
        raise ValueError("workers must be 1 or 2")
    manifest_path = manifest_path.resolve()
    output = output.resolve()
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads(manifest_path.read_text())
    jobs = _validate_manifest(manifest)
    for job in jobs:
        job_output = Path(_flag(job["argv"], "--out_folder", job["job_id"])).resolve()
        if (
            output == job_output
            or output in job_output.parents
            or job_output in output.parents
        ):
            raise ValueError(
                "Results directory overlaps a LigandMPNN job output folder"
            )
    output.mkdir(parents=True)
    statuses: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_run_job, job, output): job for job in jobs}
        for future in as_completed(futures):
            status = future.result()
            statuses[status["job_id"]] = status
    samples: list[dict[str, Any]] = []
    for job in jobs:
        for sample in audit_job(job["audit_manifest"], job):
            sample["parent"] = str(job["parent"])
            sample["arm"] = job["arm"]
            samples.append(sample)
    if len(samples) != EXPECTED_SAMPLES:
        raise ValueError(f"Auditor returned {len(samples)} samples; expected 48")
    groups: dict[str, list[str]] = {}
    by_arm: dict[str, int] = {}
    for sample in samples:
        key = sample["alpha_sequence"] + ":" + sample["beta_sequence"]
        groups.setdefault(key, []).append(sample["sample_id"])
        by_arm[sample["arm"]] = by_arm.get(sample["arm"], 0) + 1
    report: dict[str, Any] = {
        "schema_version": 1,
        "manifest": str(manifest_path),
        "source_sha256": manifest["source_sha256"],
        "jobs": [statuses[job["job_id"]] for job in jobs],
        "samples": samples,
        "duplicate_groups": [
            {"sequence": sequence, "sample_ids": ids}
            for sequence, ids in groups.items()
            if len(ids) > 1
        ],
        "summary": {
            "job_count": len(jobs),
            "raw_sequence_count": len(samples),
            "unique_sequence_count": len(groups),
            "by_arm": by_arm,
        },
    }
    _export(output, report)
    LOGGER.info("Completed %d audited designs across %d jobs", len(samples), len(jobs))
    return report
