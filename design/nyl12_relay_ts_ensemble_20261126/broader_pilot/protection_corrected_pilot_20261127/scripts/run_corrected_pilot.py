#!/usr/bin/env python3
"""Run the frozen corrected pilot only after root WT and resource releases."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
from typing import Any

PILOT = Path(__file__).resolve().parents[1]
RELEASE = PILOT / "root_wt_acceptance_release_v3.json"
MAX_PREFLIGHT_AGE_SECONDS = 1800


def digest(path: Path) -> str:
    """Return the SHA256 digest of exact file bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_release(manifest_path: Path, release_path: Path = RELEASE) -> dict[str, Any]:
    """Validate the exact WT receipt, manifest pin, and shared-host authorization."""
    release = json.loads(release_path.read_text(encoding="utf-8"))
    required = {
        "approved": True,
        "accepted_WT_source_id": "WT_paired_seed_20261126",
        "force_converged": True,
        "positive_phi_A_B_sites": 4,
        "manifest_sha256": digest(manifest_path),
        "allow_shared_host": True,
        "max_cpu_threads": 2,
    }
    for key, value in required.items():
        if release.get(key) != value:
            raise ValueError(f"Root release mismatch: {key}")
    if not release.get("gpu_uuid", "").startswith("GPU-"):
        raise ValueError("Root release must pin an explicit CUDA GPU UUID")
    receipt_path = Path(release["single_control_acceptance_path"])
    if digest(receipt_path) != release["single_control_acceptance_sha256"]:
        raise ValueError("Pinned WT acceptance receipt hash changed")
    if json.loads(receipt_path.read_text(encoding="utf-8")).get("accepted") is not True:
        raise ValueError("Pinned WT acceptance receipt is not accepted")
    return release


def verify_resource(release: dict[str, Any]) -> dict[str, Any]:
    """Require the pinned fresh Violet receipt and an idle GPU at launch time."""
    preflight_path = Path(release["resource_release_path"])
    if digest(preflight_path) != release["resource_release_sha256"]:
        raise ValueError("Pinned root resource release hash changed")
    resource = json.loads(preflight_path.read_text(encoding="utf-8"))
    required = {
        "approved": True,
        "host": "violet",
        "gpu_idle": True,
        "no_unrelated_gpu_processes": True,
        "allow_shared_host": True,
        "max_cpu_threads": 2,
        "reviewed_sessions": True,
        "reviewed_cpu_gpu_processes": True,
        "reviewed_queued_work": True,
        "scientific_jobs_must_continue": True,
        "gpu_uuid": release["gpu_uuid"],
    }
    for key, value in required.items():
        if resource.get(key) != value:
            raise ValueError(f"Resource release mismatch: {key}")
    source_path = Path(resource["source_preflight_path"])
    if digest(source_path) != resource["source_preflight_sha256"]:
        raise ValueError("Root source resource preflight hash changed")
    checked = datetime.fromisoformat(resource["created_utc"])
    age = (datetime.now(timezone.utc) - checked).total_seconds()
    if age < 0 or age > MAX_PREFLIGHT_AGE_SECONDS:
        raise ValueError("Violet resource preflight is stale")
    source = json.loads(source_path.read_text(encoding="utf-8"))
    rows = source["gpu"].splitlines()[1:]
    gpu0 = next((row for row in rows if row.split(",", 1)[0].strip() == "0"), None)
    if gpu0 is None:
        raise ValueError("Violet resource receipt has no GPU0 row")
    fields = [field.strip() for field in gpu0.split(",")]
    if len(fields) < 6 or fields[1] != release["gpu_uuid"]:
        raise ValueError("Root GPU UUID differs from the Violet GPU0 receipt")
    memory_mib = int(re.search(r"\d+", fields[3]).group())
    utilization = int(re.search(r"\d+", fields[5]).group())
    if memory_mib > 256 or utilization != 0:
        raise RuntimeError("Violet GPU0 is occupied")
    process_rows = [
        row for row in source["gpu_compute_processes"].splitlines()[1:] if row.strip()
    ]
    if process_rows:
        raise RuntimeError("Source resource receipt reports active compute processes")
    live_devices = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=uuid,memory.used,utilization.gpu",
            "--format=csv,noheader",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    device = next(
        (
            row
            for row in live_devices
            if row.split(",", 1)[0].strip() == release["gpu_uuid"]
        ),
        None,
    )
    if device is None:
        raise RuntimeError("Pinned GPU UUID is unavailable on Violet")
    live_fields = [field.strip() for field in device.split(",")]
    live_memory = int(re.search(r"\d+", live_fields[1]).group())
    live_utilization = int(re.search(r"\d+", live_fields[2]).group())
    if live_memory > 256 or live_utilization != 0:
        raise RuntimeError("Pinned Violet GPU UUID is occupied")
    active_gpu_uuids = subprocess.run(
        ["nvidia-smi", "--query-compute-apps=gpu_uuid", "--format=csv,noheader"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    if release["gpu_uuid"] in active_gpu_uuids:
        raise RuntimeError("Pinned Violet GPU UUID has an active compute process")
    return resource


def namespace_ids(
    samples: list[dict[str, Any]], candidates: list[dict[str, Any]]
) -> None:
    """Give raw samples, candidate IDs, and all cross-references a new namespace."""
    for sample in samples:
        sample["sample_id"] = f"corr01_{sample['sample_id']}"
    for candidate in candidates:
        suffix = candidate["candidate_id"].rsplit("_", 1)[-1]
        candidate["candidate_id"] = f"relay_corr01_{suffix}"
        candidate["sample_id"] = f"corr01_{candidate['sample_id']}"
        candidate["lineage"] = [f"corr01_{item}" for item in candidate["lineage"]]


def write_effective_manifest(
    manifest: dict[str, Any],
    output_path: Path,
    manifest_sha: str,
    release_path: Path,
    resource_path: Path,
    gpu_uuid: str,
    max_cpu_threads: int,
) -> dict[str, Any]:
    """Persist the resolved runtime manifest and exact release provenance."""
    with output_path.open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    effective_sha = digest(output_path)
    acceptance = {
        "manifest_sha256": manifest_sha,
        "root_wt_release_sha256": digest(release_path),
        "resource_release_sha256": digest(resource_path),
        "effective_manifest_sha256": effective_sha,
        "effective_cuda_visible_devices": gpu_uuid,
        "max_cpu_threads": max_cpu_threads,
    }
    acceptance_path = output_path.with_name(f"{output_path.stem}_acceptance.json")
    with acceptance_path.open("x", encoding="utf-8") as handle:
        json.dump(acceptance, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return acceptance


def main() -> None:
    """Run one immutable attempt after exact scientific and resource gates pass."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    release = verify_release(args.manifest)
    resource = verify_resource(release)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if (
        manifest["request"]["campaign"]
        != "nyl12_relay_protection_corrected_broad_pilot_20261127"
    ):
        raise ValueError("Unexpected campaign identity")
    environment = manifest["native"]["environment"]
    environment["CUDA_VISIBLE_DEVICES"] = release["gpu_uuid"]
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        if int(environment[variable]) > resource["max_cpu_threads"]:
            raise ValueError(f"Thread setting exceeds root release: {variable}")
    runner_path = Path(manifest["native"]["runner"])
    import sys

    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("corrected_broad_runner", runner_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load frozen source runner")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    native_pool = runner.pool_samples

    def corrected_pool(raw_samples: list[dict[str, Any]]) -> list[dict[str, Any]]:
        candidates = native_pool(raw_samples)
        namespace_ids(raw_samples, candidates)
        return candidates

    runner.pool_samples = corrected_pool
    write_effective_manifest(
        manifest,
        PILOT / f"{args.output.name}_effective_manifest.json",
        digest(args.manifest),
        RELEASE,
        Path(release["resource_release_path"]),
        release["gpu_uuid"],
        resource["max_cpu_threads"],
    )
    runner.run(manifest, args.output)


if __name__ == "__main__":
    main()
