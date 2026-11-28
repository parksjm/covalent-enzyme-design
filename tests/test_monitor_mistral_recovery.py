"""Focused checks for read-only Mistral recovery progress monitoring."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import monitor_mistral_recovery as watcher
from scripts.monitor_mistral_recovery import ATTEMPT, inspect, persist


def setup_run(root: Path, *, exit_status: int | None = None) -> Path:
    """Create one controller with a two-case benchmark and receipt index."""
    attempt = root / ATTEMPT
    (attempt / "controller_job").mkdir(parents=True)
    (attempt / "controller_job" / "status.json").write_text(
        json.dumps(
            {
                "started_utc": "2026-11-26T02:00:00+00:00",
                "exit_status": exit_status,
            }
        )
    )
    (root / "cases.json").write_text(
        json.dumps({"cases": [{"case_id": "a"}, {"case_id": "b"}]})
    )
    return root


def test_ignores_old_receipts_and_counts_new_success_and_failure(
    tmp_path: Path,
) -> None:
    """Old attempts are excluded and only complete stop receipts count valid."""
    root = setup_run(tmp_path / "run")
    model = root / "results" / "mistral"
    model.mkdir(parents=True)
    refs = {}
    for case, started, status, finish in [
        ("a", "2026-11-26T01:00:00+00:00", "complete", "stop"),
        ("b", "2026-11-26T02:01:00+00:00", "failed", None),
    ]:
        rel = f"attempts/{case}/receipt.json"
        path = model / rel
        path.parent.mkdir(parents=True)
        path.write_text(
            json.dumps(
                {
                    "started_utc": started,
                    "ended_utc": started,
                    "process_status": status,
                    "finish_reason": finish,
                }
            )
        )
        refs[case] = {"receipt": rel}
    (model / "progress.json").write_text(json.dumps({"cases": refs}))
    snapshot = inspect(root, tmp_path / "out")
    assert (snapshot["valid"], snapshot["new_failed"], snapshot["phase"]) == (
        0,
        1,
        "loading_or_warmup",
    )


def test_terminal_success_requires_all_twenty_cases(tmp_path: Path) -> None:
    """A zero controller exit with incomplete receipts is a failure."""
    root = setup_run(tmp_path / "run", exit_status=0)
    snapshot = inspect(root, tmp_path / "out")
    assert snapshot["phase"] == "failed"


def test_resume_directory_changes_phase_and_persist_deduplicates_events(
    tmp_path: Path,
) -> None:
    """Phase transition is recorded, while repeated identical states are not."""
    root = setup_run(tmp_path / "run")
    attempt = root / ATTEMPT
    (attempt / "benchmark_resume").mkdir()
    snapshot = inspect(root, tmp_path / "out")
    assert snapshot["phase"] == "canary/inference"
    persist(tmp_path / "out", snapshot)
    persist(tmp_path / "out", snapshot)
    assert len((tmp_path / "out" / "events.jsonl").read_text().splitlines()) == 1


@pytest.mark.parametrize(
    "terminal,expected", [("complete", 0), ("failed", 1), ("timeout", 2)]
)
def test_monitor_terminates_and_preserves_timeout_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, terminal: str, expected: int
) -> None:
    """Stop on terminal results or the deadline instead of polling forever."""
    phases = iter(["canary/inference", terminal])
    times = iter([0, 4000 if terminal == "timeout" else 1])
    recorded = []

    def probe(_root: Path, _output: Path) -> dict:
        return {"phase": next(phases), "valid": 0, "total": 20, "new_failed": 0}

    monkeypatch.setattr(watcher, "inspect", probe)
    monkeypatch.setattr(
        watcher, "persist", lambda _path, snapshot: recorded.append(dict(snapshot))
    )
    monkeypatch.setattr(watcher.time, "monotonic", lambda: next(times))
    monkeypatch.setattr(watcher.time, "sleep", lambda _seconds: None)
    assert watcher.monitor(tmp_path, tmp_path / "out", 60, 1) == expected
    assert recorded[-1]["phase"] == terminal
    if terminal == "timeout":
        assert recorded[-1]["timeout"] is True
