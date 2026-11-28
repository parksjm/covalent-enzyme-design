"""Completion evidence, source freezing and one-shot workflow invocation."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from nyl12_covalent import matched_wt_analysis_queue as queue
from nyl12_covalent.matched_wt_gate import ReceiptReader, digest, inspect_cases


def write(path: str | Path, record: dict[str, Any]) -> None:
    Path(path).write_text(json.dumps(record) + "\n")


@pytest.fixture
def config(tmp_path: Path) -> dict[str, Any]:
    source = tmp_path / "source.py"
    source.write_text("# frozen\n")
    manifest = tmp_path / "manifest.json"
    write(manifest, {"cases": []})
    worker = tmp_path / "worker.json"
    write(worker, {"status": "running"})
    cases = []
    for i in range(12):
        folder = tmp_path / str(i)
        folder.mkdir()
        topology = folder / "system.prmtop"
        topology.write_text("topology")
        case = {
            "id": str(i),
            "acceptance_05": str(folder / "05.json"),
            "acceptance_06": str(folder / "06.json"),
            "replica_complete": str(folder / "complete.json"),
            "replica_job_status": str(folder / "job.json"),
            "production_topology": str(topology),
            "expected_topology_sha256": digest(str(topology)),
        }
        for stage, steps, time in ((5, 5_000_000, 10500), (6, 45_000_000, 100500)):
            write(
                case[f"acceptance_0{stage}"],
                {
                    "completed_steps": steps,
                    "final_time_ps": time,
                    "final_temperature_K": 300,
                    "temperature_gate_pass": True,
                    "total_production_ns": 100,
                },
            )
        write(case["replica_complete"], {"completed": True})
        write(
            case["replica_job_status"],
            {
                "ended_utc": "2026-11-21T00:00:00Z",
                "exit_status": 0,
            },
        )
        cases.append(case)
    return {
        "cases": cases,
        "worker_statuses": [str(worker)],
        "source_sha256": {str(source): digest(str(source))},
        "source_manifest": str(manifest),
        "source_manifest_sha256": digest(str(manifest)),
        "poll_seconds": 60,
        "deadline_hours": 48,
        "command": ["explicit", "workflow"],
        "cwd": str(tmp_path),
    }


def test_runs_once_after_all_complete(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = []
    monkeypatch.setattr(
        queue.subprocess,
        "run",
        lambda *a, **k: calls.append((a, k)) or SimpleNamespace(returncode=0),
    )
    output = tmp_path / "controller"
    queue.run_queue(config, output)
    assert len(calls) == 1
    assert calls[0] == ((config["command"],), {"cwd": config["cwd"], "check": False})
    assert json.loads((output / "status.json").read_text())["status"] == "completed"
    with pytest.raises(FileExistsError):
        queue.run_queue(config, output)
    assert len(calls) == 1


def test_10ns_only_waits_without_launch(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = config["cases"][-1]
    record = Path(case["acceptance_06"]).read_text()
    Path(case["acceptance_06"]).unlink()
    sleeps = []
    calls = []

    def finish(seconds: float) -> None:
        assert not calls
        sleeps.append(seconds)
        Path(case["acceptance_06"]).write_text(record)

    monkeypatch.setattr(queue.time, "sleep", finish)
    monkeypatch.setattr(
        queue.subprocess,
        "run",
        lambda *a, **k: calls.append(a) or SimpleNamespace(returncode=0),
    )
    queue.run_queue(config, tmp_path / "controller")
    assert sleeps == [60]
    assert len(calls) == 1


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("completed_steps", 5_000_000),
        ("final_time_ps", 10500),
        ("final_time_ps", 100500.01),
        ("temperature_gate_pass", False),
        ("final_temperature_K", 351),
        ("total_production_ns", 10),
    ],
)
def test_bad_stage_rejected(config: dict, field: str, value: Any) -> None:
    path = config["cases"][0]["acceptance_06"]
    record = json.loads(Path(path).read_text())
    record[field] = value
    write(path, record)
    with pytest.raises(ValueError):
        inspect_cases(config, ReceiptReader())


def test_missing_flag_is_explicit_error(config: dict) -> None:
    path = config["cases"][0]["acceptance_05"]
    record = json.loads(Path(path).read_text())
    del record["temperature_gate_pass"]
    write(path, record)
    with pytest.raises(KeyError):
        inspect_cases(config, ReceiptReader())


@pytest.mark.parametrize(
    "kind", ["worker", "job", "topology", "source", "manifest", "workflow"]
)
def test_failure_status_no_unintended_launch(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    case = config["cases"][0]
    if kind == "worker":
        write(config["worker_statuses"][0], {"status": "failed"})
    elif kind == "job":
        write(
            case["replica_job_status"],
            {"ended_utc": "2026-11-21T00:00:00Z", "exit_status": 1},
        )
    elif kind == "topology":
        Path(case["production_topology"]).write_text("wrong")
    elif kind == "source":
        Path(next(iter(config["source_sha256"]))).write_text("wrong")
    elif kind == "manifest":
        Path(config["source_manifest"]).write_text("wrong")
    calls = []
    monkeypatch.setattr(
        queue.subprocess,
        "run",
        lambda *a, **k: calls.append(a) or SimpleNamespace(returncode=3),
    )
    output = tmp_path / "controller"
    with pytest.raises((RuntimeError, ValueError)):
        queue.run_queue(config, output)
    assert len(calls) == (1 if kind == "workflow" else 0)
    status = json.loads((output / "status.json").read_text())
    assert status["status"] == "failed"
    assert status["recovery_hint"]


def test_deadline(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    clock = iter([0, 48 * 3600])
    monkeypatch.setattr(queue.time, "monotonic", lambda: next(clock))
    with pytest.raises(TimeoutError):
        queue.run_queue(config, tmp_path / "controller")


def test_sources_verified_again_before_command(
    config: dict, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = queue.verify_topologies

    def mutate(conf: dict) -> None:
        original(conf)
        Path(conf["source_manifest"]).write_text("changed while waiting")

    monkeypatch.setattr(queue, "verify_topologies", mutate)
    with pytest.raises(ValueError, match="Frozen source"):
        queue.run_queue(config, tmp_path / "controller")


def test_receipts_pending_and_malformed(config: dict) -> None:
    case = config["cases"][0]
    Path(case["replica_complete"]).unlink()
    write(case["replica_job_status"], {"ended_utc": None, "exit_status": None})
    assert len(inspect_cases(config, ReceiptReader())) == 11
    Path(case["replica_job_status"]).write_text("{")
    reader = ReceiptReader(partial_limit=1)
    assert len(inspect_cases(config, reader)) == 11
    with pytest.raises(ValueError, match="unterminated"):
        inspect_cases(config, reader)
    Path(case["replica_job_status"]).write_text("{\n")
    with pytest.raises(json.JSONDecodeError):
        inspect_cases(config, ReceiptReader())


@pytest.mark.parametrize("record", [{"completed": False}, {}, {"completed": "true"}])
def test_terminal_completed_flag_required(config: dict, record: dict) -> None:
    write(config["cases"][0]["replica_complete"], record)
    with pytest.raises(ValueError):
        inspect_cases(config, ReceiptReader())


@pytest.mark.parametrize(
    "update",
    [
        {"cases": []},
        {"worker_statuses": []},
        {"poll_seconds": 1},
        {"deadline_hours": 0},
        {"deadline_hours": float("nan")},
        {"poll_seconds": float("inf")},
        {"command": []},
        {"command": [1]},
        {"cwd": "/does/not/exist"},
        {"source_sha256": {}},
    ],
)
def test_config_invalid(config: dict, tmp_path: Path, update: dict) -> None:
    config.update(update)
    with pytest.raises(ValueError):
        queue.run_queue(config, tmp_path / "controller")
