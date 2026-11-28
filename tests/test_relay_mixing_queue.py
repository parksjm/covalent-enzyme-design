"""Deferred six-check queue gates and process receipt isolation."""

from __future__ import annotations
import json
from pathlib import Path
from unittest.mock import Mock
import pytest
from nyl12_covalent import relay_mixing_queue as queue
from nyl12_covalent.neb_umbrella_campaign import digest, write_json


def accepted_main(tmp_path: Path) -> Path:
    """Create all57 independently accepted terminal native checkpoint receipts."""
    for label in queue.MAIN_LABELS:
        p = tmp_path / label / "14_production_5ps"
        p.mkdir(parents=True)
        (p / "restart.ncrst").write_bytes(b"accepted")
        write_json(
            p / "acceptance.json",
            dict(
                accepted=True,
                completed_steps=10000,
                every_integration_step_SCC_converged=True,
                frozen_outer_max_displacement_A=0,
                checkpoint_sha256=digest(p / "restart.ncrst"),
            ),
        )
        write_json(p / "process.json", dict(exit_status=0, ended_utc="2026-11-27"))
    return tmp_path


def test_receipt_gate(tmp_path: Path) -> None:
    """Missing or failed receipt, changed checkpoint and active worker block launch."""
    assert not queue.main_gate(tmp_path, [])["ready"]
    base = accepted_main(tmp_path)
    assert queue.main_gate(base, [])["ready"]
    assert not queue.main_gate(base, [123])["ready"]
    f = base / queue.MAIN_LABELS[0] / "14_production_5ps/process.json"
    write_json(f, dict(exit_status=1, ended_utc="now"))
    assert not queue.main_gate(base, [])["ready"]
    write_json(f, dict(exit_status=0, ended_utc="now"))
    f.with_name("restart.ncrst").write_bytes(b"changed")
    assert not queue.main_gate(base, [])["ready"]


def method() -> dict:
    """Return the sole authorized Hamiltonian inventory."""
    return dict(
        system=dict(
            full_atoms=104922,
            active_atoms=4165,
            qm_real_atoms=98,
            qm_mm_link_atoms=8,
            qm_charge=-2,
            qm_method="DFTB3/3ob-3-1",
        ),
        temperature_K=300,
        production_horizon_ps=30,
        parent_prepared="/tmp/main/prepared.json",
    )


def test_method() -> None:
    """Changed method/thermal protocol is never substituted."""
    cfg = method()
    queue.validate_method(cfg)
    cfg["temperature_K"] = 310
    with pytest.raises(ValueError):
        queue.validate_method(cfg)
    cfg["system"]["qm_charge"] = 0
    with pytest.raises(ValueError):
        queue.validate_method(cfg)


def test_six_queue(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Every unique label gets a queued status; one false poll precedes execution."""
    tasks = [dict(label=label) for label in queue.TARGETS]
    path = tmp_path / "tasks.json"
    write_json(
        path, dict(tasks=tasks, prepared="/tmp/frozen.json", prepared_sha256="x")
    )
    monkeypatch.setattr(queue, "verify_target", Mock(return_value=method()))
    monkeypatch.setattr(queue, "main_worker_pids", Mock(return_value=[]))
    monkeypatch.setattr(
        queue,
        "main_gate",
        Mock(
            side_effect=[
                dict(ready=False, accepted60ps_count=56, active_main_worker_pids=[]),
                dict(ready=True, accepted60ps_count=57, active_main_worker_pids=[]),
            ]
        ),
    )
    monkeypatch.setattr(queue.time, "sleep", Mock())
    execute = Mock()
    monkeypatch.setattr(queue, "execute_task", execute)
    out = tmp_path / "queue"
    queue.run_queue(path, out, 0)
    assert execute.call_count == 6
    assert all(
        json.loads((out / f"{label}.json").read_text())["status"] == "queued"
        for label in queue.TARGETS
    )
    assert (out / "queue_completion.json").exists()
    tasks.pop()
    write_json(path, dict(tasks=tasks))
    with pytest.raises(ValueError):
        queue.run_queue(path, tmp_path / "bad", 0)


def test_execute_task(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Per-label stdout and process failure remain independently inspectable."""
    monkeypatch.setattr(queue.subprocess, "run", Mock(return_value=Mock(returncode=0)))
    task = dict(label="A026_late", argv=["fake"], cwd=str(tmp_path), environment={})
    queue.execute_task(task, tmp_path)
    assert json.loads((tmp_path / "A026_late.json").read_text())["exit_status"] == 0
    task["label"] = "A028_late"
    monkeypatch.setattr(queue.subprocess, "run", Mock(return_value=Mock(returncode=2)))
    with pytest.raises(RuntimeError):
        queue.execute_task(task, tmp_path)
    assert json.loads((tmp_path / "A028_late.json").read_text())["status"] == "failed"


def test_cli(monkeypatch: pytest.MonkeyPatch) -> None:
    """CLI delegates its exact task and fresh output paths."""
    import sys

    runner = Mock()
    monkeypatch.setattr(queue, "run_queue", runner)
    monkeypatch.setattr(
        sys, "argv", ["queue", "--tasks", "tasks.json", "--output", "fresh"]
    )
    queue.main()
    runner.assert_called_once_with(Path("tasks.json"), Path("fresh"))
