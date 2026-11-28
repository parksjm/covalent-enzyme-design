"""Validate process identity checks and whole-campaign restoration ordering."""

from __future__ import annotations

import importlib.util
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
from types import ModuleType
from typing import Any

import pytest


@pytest.fixture
def module() -> ModuleType:
    """Load the standalone controller without launching it."""
    path = Path(__file__).parents[1] / "scripts/relay_campaign_resume.py"
    spec = importlib.util.spec_from_file_location("relay_campaign_resume", path)
    assert spec and spec.loader
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@pytest.fixture
def row() -> dict[str, Any]:
    """Provide one exact process identity."""
    return dict(pid=10, start_ticks=20, uid=30, argv="sander ", cwd="/work", state="T")


def test_matches_all_identity_fields(module: ModuleType, row: dict[str, Any]) -> None:
    """Reject PID reuse and changes in command, owner or working directory."""
    assert module.matches(row, row)
    for key in ("pid", "start_ticks", "uid", "argv", "cwd"):
        assert not module.matches(dict(row, **{key: "different"}), row)


def test_resume_validates_every_process_before_any_signal(
    module: ModuleType, row: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A mismatched second process prevents even the first SIGCONT."""
    boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    signals = []
    closed = []
    monkeypatch.setattr(module.os, "getuid", lambda: 30)
    monkeypatch.setattr(module.os, "pidfd_open", lambda pid: pid + 100, raising=False)
    monkeypatch.setattr(module.os, "close", closed.append)
    monkeypatch.setattr(
        module.signal,
        "pidfd_send_signal",
        lambda *args: signals.append(args),
        raising=False,
    )
    monkeypatch.setattr(module, "identity", lambda pid: row)
    with pytest.raises(ValueError, match="identity changed"):
        module.resume_all([row, dict(row, pid=11)], boot)
    assert not signals
    assert closed == [110, 111]


def test_resume_all_exact_processes(
    module: ModuleType, row: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Resume every validated process through its own pidfd."""
    boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    signals = []
    monkeypatch.setattr(module.os, "getuid", lambda: 30)
    monkeypatch.setattr(module.os, "pidfd_open", lambda pid: pid + 100, raising=False)
    monkeypatch.setattr(module.os, "close", lambda fd: None)
    monkeypatch.setattr(
        module.signal,
        "pidfd_send_signal",
        lambda *args: signals.append(args),
        raising=False,
    )
    monkeypatch.setattr(module, "identity", lambda pid: dict(row, pid=pid))
    assert module.resume_all([row, dict(row, pid=11)], boot) == [10, 11]
    assert signals == [(110, signal.SIGCONT), (111, signal.SIGCONT)]


def test_live_rejects_reused_or_zombie_process(
    module: ModuleType, row: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Never wait for an unrelated process that reused the wrapper PID."""
    monkeypatch.setattr(module, "identity", lambda pid: row)
    assert module.live(row)
    assert not module.live(dict(row, start_ticks=21))
    monkeypatch.setattr(module, "identity", lambda pid: dict(row, state="Z"))
    assert not module.live(row)


def test_run_restores_once_only_after_all_checks_idle(
    module: ModuleType,
    row: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The controller does not restore a lone early-finishing worker."""
    events = []
    states = iter([True, False, False, False])
    monkeypatch.setattr(module.os, "pidfd_open", lambda pid: pid, raising=False)
    monkeypatch.setattr(
        module.signal, "pidfd_send_signal", lambda *args: None, raising=False
    )
    monkeypatch.setattr(module, "live", lambda _: next(states))
    monkeypatch.setattr(module, "target_busy", lambda _: False)
    monkeypatch.setattr(module, "journal", lambda *args: events.append(args[1]))
    monkeypatch.setattr(module.time, "sleep", lambda _: events.append("wait"))
    monkeypatch.setattr(
        module, "resume_all", lambda *_: events.append("resume") or [10]
    )
    config = dict(
        status=str(tmp_path / "status.json"),
        targets=[
            dict(label=str(i), root=str(tmp_path / str(i)), wrapper=row)
            for i in range(2)
        ],
        poll_seconds=60,
        processes=[row],
        boot_id="boot",
        known_wrapper_issue="retained legacy wrapper",
    )
    module.run(config)
    assert events == ["intent", "wait", "intent", "resume", "outcome"]


def test_real_identity_and_atomic_save(module: ModuleType, tmp_path: Path) -> None:
    """Read this test's process and atomically preserve an operational receipt."""
    row = module.identity(os.getpid())
    assert row["uid"] == os.getuid()
    assert row["start_ticks"] > 0
    path = tmp_path / "state.json"
    module.save(path, row)
    assert json.loads(path.read_text()) == row
    assert not path.with_suffix(".tmp").exists()
    assert not module.target_busy(tmp_path)


def test_live_absent_process(
    module: ModuleType, row: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """An already-exited wrapper does not prevent campaign restoration."""

    def absent(pid: int) -> None:
        raise FileNotFoundError(pid)

    monkeypatch.setattr(module, "identity", absent)
    assert not module.live(row)


def test_boot_and_owner_guards(
    module: ModuleType, row: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject changed boot and foreign owner before sending signals."""
    with pytest.raises(ValueError, match="Boot changed"):
        module.resume_all([row], "wrong")
    monkeypatch.setattr(module.os, "getuid", lambda: 31)
    boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    with pytest.raises(ValueError, match="different owner"):
        module.resume_all([row], boot)


def test_main_hash_guard_and_dispatch(
    module: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Only a pinned control manifest reaches the monitor."""
    path = tmp_path / "config.json"
    path.write_text('{"test": true}')
    args = ["resume", "--config", str(path), "--sha256", "wrong"]
    monkeypatch.setattr(sys, "argv", args)
    with pytest.raises(ValueError, match="hash changed"):
        module.main()
    args[-1] = hashlib.sha256(path.read_bytes()).hexdigest()
    calls = []
    monkeypatch.setattr(module, "run", calls.append)
    module.main()
    assert calls == [{"test": True}]


def test_journal_requires_success(
    module: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Audit recording uses an explicit system interpreter and checks exit status."""
    calls = []
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: calls.append((a, k)))
    module.journal(dict(logger="logger.py", journal="actions.jsonl"), "intent", "test")
    assert calls[0][0][0][:2] == ["/usr/bin/python3", "logger.py"]
    assert calls[0][1]["check"] is True


def test_run_rejects_unsupported_interpreter(
    module: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Missing process-control APIs fail before monitoring or restoration."""
    monkeypatch.delattr(module.os, "pidfd_open", raising=False)
    with pytest.raises(RuntimeError, match="system Python"):
        module.run({})
