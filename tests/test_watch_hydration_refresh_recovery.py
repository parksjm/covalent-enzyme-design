"""Tests for fail-closed hydration postprocess recovery."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from scripts import watch_hydration_refresh_recovery as watcher


def put_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def fixture_block(tmp_path: Path) -> tuple[dict, dict, Path]:
    watcher.BASE = tmp_path / "campaign"
    block = watcher.BASE / "extend_02"
    prior = watcher.BASE / "extend_01" / "restart.ncrst"
    prior.parent.mkdir(parents=True)
    prior.write_bytes(b"predecessor restart")
    start_hash = hashlib.sha256(prior.read_bytes()).hexdigest()
    put_json(
        block / "status.json",
        {
            "status": "accepted",
            "exit_code": 0,
            "predecessor": "extend_01",
            "start_restart_sha256": start_hash,
        },
    )
    put_json(
        block / "process.json",
        {
            "exit_code": 0,
            "start_restart_sha256": start_hash,
        },
    )
    put_json(
        block / "scientific_acceptance.json",
        {
            "scc_verified": True,
            "native_process_passed": True,
            "finite_temperature_passed": True,
            "acyl_retained": True,
            "boundary": {"box_unchanged": True, "frozen_max_displacement_A": 5e-7},
        },
    )
    source = tmp_path / "pinned.py"
    source.write_text("pinned source\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    snapshot = watcher.BASE / "extension_source_snapshot" / "pinned.py"
    snapshot.parent.mkdir(parents=True)
    snapshot.write_text("snapshot\n")
    snap_digest = hashlib.sha256(snapshot.read_bytes()).hexdigest()
    hook = {
        "status": "failed",
        "exit_code": 1,
        "argv": ["original-hook"],
        "source_sha256_before": {str(source): digest},
        "source_sha256_after": {str(source): digest},
    }
    put_json(block / "postprocess.json", hook)
    config = {
        "segment_names": ["extend_01", "extend_02"],
        "source_snapshot_sha256": {"pinned.py": snap_digest},
    }
    recovery_source = tmp_path / "recovery.py"
    recovery_source.write_text("recovery source\n")
    recovery = {
        "id": "schema-fixed-copy",
        "command": ["python", str(recovery_source)],
        "environment_overrides": {"PYTHONPATH": "/tmp/src:/tmp"},
        "source_sha256": {
            str(recovery_source): hashlib.sha256(
                recovery_source.read_bytes()
            ).hexdigest()
        },
    }
    return config, recovery, block


def test_hash_verification_fails_on_source_change(tmp_path: Path) -> None:
    source = tmp_path / "source.py"
    source.write_text("before")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    watcher.verify_hashes({str(source): digest})
    source.write_text("after")
    with pytest.raises(ValueError, match="Recovery source changed"):
        watcher.verify_hashes({str(source): digest})


def test_eligible_waits_for_native_acceptance(tmp_path: Path) -> None:
    config, _, block = fixture_block(tmp_path)
    status = json.loads((block / "status.json").read_text())
    status["status"] = "running"
    put_json(block / "status.json", status)
    assert watcher.eligible(2, config) is None


def test_eligible_waits_for_missing_and_running_postprocess_receipts(
    tmp_path: Path,
) -> None:
    config, _, block = fixture_block(tmp_path)
    (block / "postprocess.json").unlink()
    assert watcher.eligible(2, config) is None
    put_json(block / "postprocess.json", {"status": "running"})
    assert watcher.eligible(2, config) is None


def test_eligible_errors_on_explicit_native_failure(tmp_path: Path) -> None:
    config, _, block = fixture_block(tmp_path)
    status = json.loads((block / "status.json").read_text())
    status["status"] = "failed"
    put_json(block / "status.json", status)
    with pytest.raises(ValueError, match="native segment ended as failed"):
        watcher.eligible(2, config)


def test_eligible_rejects_failed_scientific_gate(tmp_path: Path) -> None:
    config, _, block = fixture_block(tmp_path)
    acceptance = json.loads((block / "scientific_acceptance.json").read_text())
    acceptance["acyl_retained"] = False
    put_json(block / "scientific_acceptance.json", acceptance)
    with pytest.raises(ValueError, match="scientific acceptance"):
        watcher.eligible(2, config)


def test_eligible_skips_original_success_and_successful_prior_recovery(
    tmp_path: Path,
) -> None:
    config, _, block = fixture_block(tmp_path)
    hook = json.loads((block / "postprocess.json").read_text())
    hook["status"] = "completed"
    hook["exit_code"] = 0
    put_json(block / "postprocess.json", hook)
    assert watcher.eligible(2, config) == (block, [])
    hook["status"] = "failed"
    hook["exit_code"] = 1
    put_json(block / "postprocess.json", hook)
    put_json(block / "postprocess_recovery_attempt2.json", {"exit_code": 0})
    assert watcher.eligible(2, config) == (block, [])


def test_source_verification_detects_changed_snapshot(tmp_path: Path) -> None:
    config, _, block = fixture_block(tmp_path)
    hook = json.loads((block / "postprocess.json").read_text())
    watcher.verify_sources(config, hook)
    (watcher.BASE / "extension_source_snapshot/pinned.py").write_text("changed")
    with pytest.raises(ValueError, match="Pinned snapshot changed"):
        watcher.verify_sources(config, hook)


def test_recover_runs_reviewed_command_and_records_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config, recovery, block = fixture_block(tmp_path)
    monkeypatch.setattr(watcher, "ROOT", tmp_path)
    monkeypatch.setattr(watcher, "journal_action", lambda *args: None)
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(watcher.subprocess, "run", run)
    assert watcher.recover(2, config, recovery)
    assert calls[0][0] == recovery["command"]
    assert calls[0][1]["env"]["PYTHONPATH"] == "/tmp/src:/tmp"
    receipt = json.loads(
        (block / "postprocess_environment_recovery_watcher.json").read_text()
    )
    assert receipt["status"] == "succeeded"
    with pytest.raises(FileExistsError):
        watcher.recover(2, config, recovery)


def test_recover_records_failure_and_stops(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config, recovery, block = fixture_block(tmp_path)
    monkeypatch.setattr(watcher, "ROOT", tmp_path)
    monkeypatch.setattr(watcher, "journal_action", lambda *args: None)
    monkeypatch.setattr(
        watcher.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=7)
    )
    with pytest.raises(watcher.subprocess.CalledProcessError):
        watcher.recover(2, config, recovery)
    receipt = json.loads(
        (block / "postprocess_environment_recovery_watcher.json").read_text()
    )
    assert receipt["status"] == "failed"


def test_watch_times_out_without_touching_native_jobs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    watcher.BASE = tmp_path / "campaign"
    watcher.BASE.mkdir()
    put_json(
        watcher.BASE / "extension_20ps_config.json",
        {"segment_names": ["extend_01", "extend_02"]},
    )
    recovery = tmp_path / "recovery.json"
    put_json(recovery, {})
    with pytest.raises(TimeoutError, match="extend_02"):
        watcher.watch(recovery, start_index=2, timeout_seconds=0)


def test_main_passes_bounded_options_to_watcher(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    received = []
    recovery = tmp_path / "recovery.json"
    monkeypatch.setattr(
        sys,
        "argv",
        ["watch", "--recovery-config", str(recovery), "--poll-seconds", "9"],
    )
    monkeypatch.setattr(watcher, "watch", lambda *args: received.append(args))
    watcher.main()
    assert received == [(recovery, 2, 9, 86400)]
