#!/usr/bin/env python3
"""Run a disjoint manifest of one-sample Nyl12 Boltz predictions."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import queue
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class Job:
    """One explicitly seeded prediction assigned to one host."""

    campaign: str
    sequence_name: str
    model_index: int
    seed: int
    yaml_path: Path


def read_jobs(path: Path) -> list[Job]:
    """Read and validate a distributed job manifest.

    Args:
        path: CSV manifest written by the partition builder.

    Returns:
        Ordered prediction jobs.
    """
    jobs: list[Job] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            job = Job(
                campaign=row["campaign"],
                sequence_name=row["sequence_name"],
                model_index=int(row["model_index"]),
                seed=int(row["seed"]),
                yaml_path=Path(row["yaml_path"]),
            )
            if not job.yaml_path.is_file():
                raise FileNotFoundError(job.yaml_path)
            jobs.append(job)
    if not jobs:
        raise ValueError(f"No jobs in {path}")
    keys = {(job.sequence_name, job.model_index) for job in jobs}
    if len(keys) != len(jobs):
        raise ValueError(f"Duplicate sequence/model jobs in {path}")
    return jobs


def canonical_dir(output_root: Path, job: Job) -> Path:
    """Return the canonical export directory for one sequence."""
    return output_root / "canonical" / job.sequence_name


def raw_prediction_dir(output_root: Path, job: Job) -> Path:
    """Return the raw Boltz prediction directory for one job."""
    return (
        output_root
        / "runs"
        / job.sequence_name
        / f"model_{job.model_index}"
        / "output"
        / f"boltz_results_{job.sequence_name}"
        / "predictions"
        / job.sequence_name
    )


def required_paths(
    base: Path, sequence_name: str, model_index: int
) -> tuple[Path, Path]:
    """Return required structure and confidence paths for a model."""
    stem = f"{sequence_name}_model_{model_index}"
    return base / f"{stem}.pdb", base / f"confidence_{stem}.json"


def output_complete(output_root: Path, job: Job) -> bool:
    """Return whether a canonical export is complete."""
    return all(
        path.is_file()
        for path in required_paths(
            canonical_dir(output_root, job), job.sequence_name, job.model_index
        )
    )


def raw_complete(output_root: Path, job: Job) -> bool:
    """Return whether raw model_0 output is ready to merge."""
    return all(
        path.is_file()
        for path in required_paths(
            raw_prediction_dir(output_root, job), job.sequence_name, 0
        )
    )


def merge_output(output_root: Path, job: Job) -> None:
    """Merge raw model_0 artifacts into the canonical export tree."""
    source_dir = raw_prediction_dir(output_root, job)
    destination = canonical_dir(output_root, job)
    destination.mkdir(parents=True, exist_ok=True)
    for source in source_dir.glob("*model_0*"):
        target_name = source.name.replace("model_0", f"model_{job.model_index}")
        shutil.copy2(source, destination / target_name)
    auxiliary_artifacts = {
        "affinity": ".json",
        "pre_affinity": ".npz",
    }
    for prefix, suffix in auxiliary_artifacts.items():
        source = source_dir / f"{prefix}_{job.sequence_name}{suffix}"
        if source.is_file():
            target = destination / (
                f"{prefix}_{job.sequence_name}_model_{job.model_index}{suffix}"
            )
            shutil.copy2(source, target)
    if not output_complete(output_root, job):
        raise FileNotFoundError(f"Canonical merge incomplete for {job}")


def command(args: argparse.Namespace, job: Job) -> list[str]:
    """Build the provenance-stable one-sample Boltz command."""
    run_output = (
        args.output_root
        / "runs"
        / job.sequence_name
        / f"model_{job.model_index}"
        / "output"
    )
    return [
        str(args.boltz_python),
        str(args.source_dir / "run_boltz_cli.py"),
        "predict",
        str(job.yaml_path),
        "--model",
        "boltz2",
        "--output_format",
        "pdb",
        "--cache",
        str(args.cache),
        "--recycling_steps",
        "3",
        "--sampling_steps",
        "200",
        "--diffusion_samples",
        "1",
        "--seed",
        str(job.seed),
        "--out_dir",
        str(run_output),
        "--use_potentials",
        "--no_kernels",
    ]


def run_job(args: argparse.Namespace, job: Job, gpu_queue: queue.Queue[str]) -> None:
    """Run or resume one manifest job on the next available GPU."""
    if output_complete(args.output_root, job):
        LOGGER.info("Skipping complete %s model_%d", job.sequence_name, job.model_index)
        return
    if raw_complete(args.output_root, job):
        merge_output(args.output_root, job)
        LOGGER.info("Merged existing %s model_%d", job.sequence_name, job.model_index)
        return
    gpu = gpu_queue.get()
    try:
        log_path = (
            args.output_root
            / "logs"
            / f"{job.sequence_name}_model_{job.model_index}.log"
        )
        log_path.parent.mkdir(parents=True, exist_ok=True)
        invocation = command(args, job)
        environment = os.environ.copy()
        environment["CUDA_VISIBLE_DEVICES"] = gpu
        with log_path.open("a", encoding="utf-8") as log:
            log.write("Command: " + " ".join(invocation) + "\n")
            log.write(f"CUDA_VISIBLE_DEVICES={gpu}\n")
            log.flush()
            subprocess.run(
                invocation,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=environment,
            )
        merge_output(args.output_root, job)
        LOGGER.info("Completed %s model_%d", job.sequence_name, job.model_index)
    finally:
        gpu_queue.put(gpu)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--gpu-devices", default="0,1")
    parser.add_argument("--boltz-python", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    """Run all disjoint jobs and write a machine-readable completion marker."""
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    jobs = read_jobs(args.jobs)
    gpu_queue: queue.Queue[str] = queue.Queue()
    for gpu in (item.strip() for item in args.gpu_devices.split(",")):
        if gpu:
            gpu_queue.put(gpu)
    if gpu_queue.empty():
        raise ValueError("At least one GPU device is required")
    args.output_root.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Starting %d disjoint jobs on %d GPUs", len(jobs), gpu_queue.qsize())
    with ThreadPoolExecutor(max_workers=gpu_queue.qsize()) as executor:
        futures = [executor.submit(run_job, args, job, gpu_queue) for job in jobs]
        for future in as_completed(futures):
            future.result()
    marker = {
        "jobs": len(jobs),
        "complete": sum(output_complete(args.output_root, job) for job in jobs),
    }
    if marker["complete"] != marker["jobs"]:
        raise RuntimeError(f"Worker completeness failure: {marker}")
    (args.output_root / "WORKER_COMPLETE.json").write_text(
        json.dumps(marker, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info("Worker complete: %s", marker)


if __name__ == "__main__":
    main()
