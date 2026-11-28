"""Mocked resource and main-flow coverage for the corrected-pilot gate."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import run_corrected_pilot as pilot  # noqa: E402


def digest(path: Path) -> str:
    """Return exact fixture file SHA256."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resource_bundle(
    tmp_path: Path, *, created_utc: datetime | None = None
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Create a root release and source resource receipt for mocked checks."""
    source = tmp_path / "source_preflight.json"
    source_data = {
        "gpu": (
            "index, uuid, name, memory.used [MiB], memory.total [MiB], "
            "utilization.gpu [%]\n0, GPU-test, NVIDIA, 38 MiB, 49140 MiB, 0 %\n"
        ),
        "gpu_compute_processes": "gpu_uuid, pid, process_name, used_gpu_memory [MiB]\n",
    }
    source.write_text(json.dumps(source_data), encoding="utf-8")
    resource = {
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
        "gpu_uuid": "GPU-test",
        "created_utc": (created_utc or datetime.now(timezone.utc)).isoformat(),
        "source_preflight_path": str(source),
        "source_preflight_sha256": digest(source),
    }
    resource_path = tmp_path / "resource.json"
    resource_path.write_text(json.dumps(resource), encoding="utf-8")
    release = {
        "gpu_uuid": "GPU-test",
        "resource_release_path": str(resource_path),
        "resource_release_sha256": digest(resource_path),
    }
    return release, resource


def fake_nvidia_smi(memory: str = "38 MiB", utilization: str = "0 %") -> Any:
    """Return subprocess outputs for idle GPU inventory and no compute apps."""

    def run(command: list[str], **_: Any) -> SimpleNamespace:
        if any(argument.startswith("--query-gpu") for argument in command):
            return SimpleNamespace(stdout=f"GPU-test, {memory}, {utilization}\n")
        return SimpleNamespace(stdout="")

    return run


def test_verify_resource_accepts_fresh_idle_gpu(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    release, _ = resource_bundle(tmp_path)
    monkeypatch.setattr(pilot.subprocess, "run", fake_nvidia_smi())

    resource = pilot.verify_resource(release)

    assert resource["gpu_uuid"] == "GPU-test"


def test_verify_resource_rejects_source_hash_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    release, resource = resource_bundle(tmp_path)
    Path(resource["source_preflight_path"]).write_text("changed", encoding="utf-8")
    monkeypatch.setattr(pilot.subprocess, "run", fake_nvidia_smi())

    with pytest.raises(ValueError, match="source resource preflight hash changed"):
        pilot.verify_resource(release)


def test_verify_resource_rejects_stale_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stale = datetime.now(timezone.utc) - timedelta(
        seconds=pilot.MAX_PREFLIGHT_AGE_SECONDS + 1
    )
    release, _ = resource_bundle(tmp_path, created_utc=stale)
    monkeypatch.setattr(pilot.subprocess, "run", fake_nvidia_smi())

    with pytest.raises(ValueError, match="resource preflight is stale"):
        pilot.verify_resource(release)


def test_verify_resource_rejects_live_gpu_busy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    release, _ = resource_bundle(tmp_path)
    monkeypatch.setattr(
        pilot.subprocess,
        "run",
        fake_nvidia_smi(memory="1024 MiB", utilization="42 %"),
    )

    with pytest.raises(RuntimeError, match="GPU UUID is occupied"):
        pilot.verify_resource(release)


def test_main_writes_effective_manifest_and_uses_resolved_uuid(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pilot_root = tmp_path / "pilot"
    pilot_root.mkdir()
    manifest_path = tmp_path / "manifest.json"
    manifest = {
        "request": {
            "campaign": "nyl12_relay_protection_corrected_broad_pilot_20261127"
        },
        "native": {
            "runner": str(tmp_path / "frozen_runner.py"),
            "environment": {
                "CUDA_VISIBLE_DEVICES": "ROOT_RESOURCE_RELEASE_REQUIRED",
                "OMP_NUM_THREADS": "2",
                "OPENBLAS_NUM_THREADS": "2",
                "MKL_NUM_THREADS": "2",
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    wt_release = tmp_path / "wt_release.json"
    wt_release.write_text("{}", encoding="utf-8")
    resource_release = tmp_path / "resource_release.json"
    resource_release.write_text("{}", encoding="utf-8")
    release = {"gpu_uuid": "GPU-test", "resource_release_path": str(resource_release)}
    resource = {"max_cpu_threads": 2}
    observed: dict[str, Any] = {}
    candidate = {
        "candidate_id": "relay_01",
        "sample_id": "parent_1",
        "lineage": ["parent_1"],
    }
    runner = SimpleNamespace(pool_samples=lambda _: [candidate])

    def run_job(runtime_manifest: dict[str, Any], output: Path) -> None:
        raw = [{"sample_id": "parent_1", "sequence": "AAA"}]
        runner.pool_samples(raw)
        observed["candidate"] = candidate
        observed["manifest"] = runtime_manifest
        observed["output"] = output

    runner.run = run_job
    loader = SimpleNamespace(exec_module=lambda _: None)
    spec = SimpleNamespace(loader=loader)
    monkeypatch.setattr(pilot, "PILOT", pilot_root)
    monkeypatch.setattr(pilot, "RELEASE", wt_release)
    monkeypatch.setattr(pilot, "verify_release", lambda *_: release)
    monkeypatch.setattr(pilot, "verify_resource", lambda *_: resource)
    monkeypatch.setattr(
        pilot.importlib.util, "spec_from_file_location", lambda *_: spec
    )
    monkeypatch.setattr(pilot.importlib.util, "module_from_spec", lambda *_: runner)
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_corrected_pilot.py", str(manifest_path), str(tmp_path / "run01")],
    )

    pilot.main()

    assert (
        observed["manifest"]["native"]["environment"]["CUDA_VISIBLE_DEVICES"]
        == "GPU-test"
    )
    assert observed["candidate"]["candidate_id"] == "relay_corr01_01"
    assert observed["output"] == tmp_path / "run01"
    effective = json.loads((pilot_root / "run01_effective_manifest.json").read_text())
    assert effective["native"]["environment"]["CUDA_VISIBLE_DEVICES"] == "GPU-test"
    assert (pilot_root / "run01_effective_manifest_acceptance.json").exists()


def test_main_rejects_thread_setting_above_root_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "request": {
                    "campaign": "nyl12_relay_protection_corrected_broad_pilot_20261127"
                },
                "native": {
                    "environment": {
                        "CUDA_VISIBLE_DEVICES": "pending",
                        "OMP_NUM_THREADS": "3",
                        "OPENBLAS_NUM_THREADS": "2",
                        "MKL_NUM_THREADS": "2",
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(pilot, "verify_release", lambda *_: {"gpu_uuid": "GPU-test"})
    monkeypatch.setattr(pilot, "verify_resource", lambda *_: {"max_cpu_threads": 2})
    monkeypatch.setattr(
        sys, "argv", ["run_corrected_pilot.py", str(manifest_path), "unused"]
    )

    with pytest.raises(
        ValueError, match="Thread setting exceeds root release: OMP_NUM_THREADS"
    ):
        pilot.main()
