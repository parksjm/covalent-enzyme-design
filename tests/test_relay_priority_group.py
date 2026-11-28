"""Protect exact group scheduling and restoration after native-task failures."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import sys
from types import ModuleType
from typing import Any

import pytest


@pytest.fixture
def module(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """Import the standalone scheduler without executing its entry point."""
    scripts = Path(__file__).parents[1] / "scripts"
    monkeypatch.syspath_prepend(str(scripts))
    spec = importlib.util.spec_from_file_location(
        "priority_group_test", scripts / "run_relay_priority_group.py"
    )
    assert spec and spec.loader
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@pytest.fixture
def config(tmp_path: Path) -> dict[str, Any]:
    """Create a fully explicit one-task plan with a 42-identity reservation."""
    pinned = tmp_path / "input.json"
    pinned.write_text("{}")
    return {
        "sha256": {str(pinned): hashlib.sha256(pinned.read_bytes()).hexdigest()},
        "processes": [
            dict(pid=i, uid=30, start_ticks=40, argv="sander", cwd="/work")
            for i in range(42)
        ],
        "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
        "status": str(tmp_path / "status.json"),
        "logger": "/logger",
        "journal": "/journal",
        "tasks": [
            {
                "label": "target",
                "argv": ["/usr/bin/python3", "worker.py"],
                "cwd": str(tmp_path),
                "log": str(tmp_path / "worker.log"),
                "receipt": str(tmp_path / "receipt.json"),
                "environment": {"OMP_NUM_THREADS": "1"},
            }
        ],
    }


def test_verify_and_preflight(
    module: ModuleType, config: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject changed files, duplicate tasks and attempted outputs before pausing."""
    monkeypatch.setattr(module.os, "pidfd_open", lambda _: 1, raising=False)
    monkeypatch.setattr(
        module.signal, "pidfd_send_signal", lambda *_: None, raising=False
    )
    module.preflight(config)
    for files in ({}, {"relative": "abc"}, {next(iter(config["sha256"])): "wrong"}):
        with pytest.raises(ValueError):
            module.verify_files(files)
    with pytest.raises(ValueError, match="uniquely"):
        module.preflight(dict(config, tasks=config["tasks"] * 2))
    with pytest.raises(ValueError, match="complete"):
        module.preflight(dict(config, processes=[]))
    task = config["tasks"][0]
    with pytest.raises(ValueError, match="absolute"):
        module.preflight(dict(config, tasks=[dict(task, argv=["python"])]))
    with pytest.raises(ValueError, match="directory"):
        module.preflight(
            dict(config, tasks=[dict(task, cwd="/missing-priority-directory")])
        )
    Path(task["log"]).touch()
    with pytest.raises(FileExistsError):
        module.preflight(config)
    monkeypatch.delattr(module.os, "pidfd_open")
    with pytest.raises(RuntimeError, match="system Python"):
        module.preflight(config)


def test_pause_identity_and_partial_rollback(
    module: ModuleType, config: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Validate all identities before signaling; restore a partially paused group."""
    rows = config["processes"][:2]
    sent: list[tuple[int, int]] = []
    closed: list[int] = []
    monkeypatch.setattr(module.os, "getuid", lambda: 30)
    monkeypatch.setattr(module.os, "pidfd_open", lambda pid: pid + 100, raising=False)
    monkeypatch.setattr(module.os, "close", closed.append)
    monkeypatch.setattr(module, "identity", lambda pid: dict(rows[pid], state="R"))

    def send(fd: int, sig: int) -> None:
        if fd == 101 and sig == signal.SIGSTOP:
            raise ProcessLookupError("exited")
        sent.append((fd, sig))

    monkeypatch.setattr(module.signal, "pidfd_send_signal", send, raising=False)
    with pytest.raises(ProcessLookupError):
        module.pause_all(rows, config["boot_id"])
    assert sent == [(100, signal.SIGSTOP), (100, signal.SIGCONT)]
    assert closed == [100, 101]
    with pytest.raises(ValueError, match="Boot"):
        module.pause_all(rows, "wrong")
    with pytest.raises(ValueError, match="owner"):
        module.pause_all([dict(rows[0], uid=31)], config["boot_id"])
    monkeypatch.setattr(
        module, "identity", lambda _: dict(rows[0], start_ticks=99, state="R")
    )
    sent.clear()
    with pytest.raises(ValueError, match="identity"):
        module.pause_all(rows, config["boot_id"])
    assert not sent
    monkeypatch.setattr(module, "identity", lambda pid: dict(rows[pid], state="R"))
    monkeypatch.setattr(
        module.signal, "pidfd_send_signal", lambda *args: sent.append(args)
    )
    assert module.pause_all(rows, config["boot_id"]) == [0, 1]


@pytest.mark.parametrize("mode", ["success", "failed", "malformed", "launch_error"])
def test_run_restores_complete_group(
    module: ModuleType,
    config: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
) -> None:
    """Restore capacity only after every successfully launched task is idle."""
    events: list[str] = []
    monkeypatch.setattr(module, "preflight", lambda _: events.append("preflight"))
    monkeypatch.setattr(module, "verify_files", lambda _: None)
    monkeypatch.setattr(module, "journal", lambda _, phase, text: events.append(phase))
    monkeypatch.setattr(module, "pause_all", lambda *_: events.append("pause") or [1])
    monkeypatch.setattr(module, "resume_all", lambda *_: events.append("resume") or [1])

    class Child:
        pid = 123

        def wait(self) -> int:
            events.append("wait")
            path = Path(config["tasks"][0]["receipt"])
            path.write_text(
                "invalid"
                if mode == "malformed"
                else json.dumps({"status": "completed"})
            )
            return 1 if mode == "failed" else 0

    def start(*args: Any, **kwargs: Any) -> Child:
        events.append("launch")
        assert kwargs["env"]["OMP_NUM_THREADS"] == "1"
        if mode == "launch_error":
            raise OSError("cannot launch")
        return Child()

    monkeypatch.setattr(module.subprocess, "Popen", start)
    if mode == "success":
        module.run(config)
    else:
        with pytest.raises((OSError, RuntimeError, ValueError)):
            module.run(config)
    result = json.loads(Path(config["status"]).read_text())
    assert result["resumed_pids"] == [1]
    assert result["status"] == ("completed" if mode == "success" else "failed")
    assert not result["sampling_converged"]
    assert events.index("pause") < events.index("launch") < events.index("resume")
    if mode != "launch_error":
        assert events.index("wait") < events.index("resume")


def test_journal_and_cli(
    module: ModuleType,
    config: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Use the mandatory journal and reject an altered top-level manifest."""
    calls: list[Any] = []
    monkeypatch.setattr(
        module.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs))
    )
    module.journal(config, "intent", "test action")
    assert calls[0][1]["check"]
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    monkeypatch.setattr(
        sys, "argv", ["controller", "--config", str(path), "--sha256", digest]
    )
    monkeypatch.setattr(module, "run", lambda cfg: calls.append(cfg))
    module.main()
    assert calls[-1] == config
    monkeypatch.setattr(
        sys, "argv", ["controller", "--config", str(path), "--sha256", "bad"]
    )
    with pytest.raises(ValueError, match="hash"):
        module.main()
