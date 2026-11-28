"""Verify receipt-last publication and explicit collector failures."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SPEC = importlib.util.spec_from_file_location(
    "boundary_return", Path(__file__).with_name("return_accepted.py")
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.fixture
def collector(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple:
    """Create a complete-stage inventory and isolated fake file transport."""
    base = tmp_path / "native"
    base.mkdir()
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.setattr(MODULE, "BASE", base)
    monkeypatch.setattr(MODULE, "WORK", work)
    actions = []
    copies = []
    monkeypatch.setattr(MODULE, "journal", lambda *args: actions.append(args))
    payload = b"q_acyl_A,q_act_A,q_del_A\n-2,-0.8,-2.3\n"
    digest = hashlib.sha256(payload).hexdigest()
    row = {
        "stage": "window/04_production_5ps",
        "marker": "window/04_production_5ps/analysis_acceptance.json",
        "files": ["window/04_production_5ps/raw_features.csv"],
        "csv_sha256": digest,
    }
    inventory = {
        "stages": [row],
        "live": ["window/window_status.json"],
        "controller": {"exit_status": None},
    }
    monkeypatch.setattr(
        MODULE.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout=json.dumps(inventory)),
    )

    def sync(files: list[str]) -> None:
        copies.append(list(files))
        stage = base / row["stage"]
        stage.mkdir(parents=True, exist_ok=True)
        if row["marker"] in files:
            assert (stage / "raw_features.csv").read_bytes() == payload
            (base / row["marker"]).write_text(
                json.dumps({"accepted": True, "csv_sha256": digest})
            )
        elif row["files"][0] in files:
            (stage / "raw_features.csv").write_bytes(payload)
            (stage / "process.json").write_text(json.dumps({"exit_status": 0}))

    monkeypatch.setattr(MODULE, "sync_files", sync)
    return base, work, actions, copies, row, inventory


def test_receipt_published_after_verified_payload(collector: tuple) -> None:
    """A stage becomes eligible only after complete payload checks."""
    base, work, actions, copies, row, _ = collector
    assert not MODULE.poll()
    assert copies[:2] == [row["files"], [row["marker"]]]
    assert (base / row["marker"]).exists()
    assert json.loads((work / "return_state.json").read_text())["new_stages"] == [
        row["stage"]
    ]
    assert actions[0][0] == "intent" and actions[-1][0] == "outcome"


def test_bad_csv_never_publishes_acceptance(collector: tuple) -> None:
    """A corrupt transferred scientific payload cannot gain an acceptance marker."""
    base, _, _, copies, row, _ = collector
    row["csv_sha256"] = "wrong"
    with pytest.raises(ValueError, match="CSV identity mismatch"):
        MODULE.poll()
    assert copies == [row["files"]]
    assert not (base / row["marker"]).exists()


def test_terminal_controller_failure_is_explicit(collector: tuple) -> None:
    """Native failures preserve copied artifacts and stop the return monitor."""
    base, _, _, _, row, inventory = collector
    inventory["controller"]["exit_status"] = 1
    with pytest.raises(RuntimeError, match="controller failed"):
        MODULE.poll()
    assert (base / row["marker"]).exists()


def test_completed_stage_is_not_retransferred(collector: tuple) -> None:
    """Repeated polls preserve immutable accepted payloads."""
    _, _, _, copies, row, inventory = collector
    assert not MODULE.poll()
    copies.clear()
    inventory["controller"]["exit_status"] = 0
    assert MODULE.poll()
    assert copies == [inventory["live"]]


def test_transport_uses_exact_files_without_delete(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Transport copies only explicit paths and never deletes remote/local data."""
    calls = []
    monkeypatch.setattr(
        MODULE.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs))
    )
    MODULE.sync_files([])
    assert calls == []
    MODULE.sync_files(["window/stage/raw_features.csv"])
    command = calls[0][0][0]
    assert "--files-from=-" in command and "--delete" not in command
    assert calls[0][1]["input"] == "window/stage/raw_features.csv\n"
    assert calls[0][1]["check"] is True


def test_journal_errors_are_not_suppressed(monkeypatch: pytest.MonkeyPatch) -> None:
    """The logging prerequisite uses checked execution and records result text."""
    calls = []
    monkeypatch.setattr(
        MODULE.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs))
    )
    MODULE.journal("outcome", "test action", "test result")
    assert calls[0][1]["check"] is True
    assert calls[0][0][0][-2:] == ["--result", "test result"]


def test_main_records_failure_and_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """An SSH/transport error cannot leave a success-shaped monitor outcome."""
    actions = []
    monkeypatch.setattr(MODULE, "journal", lambda *args: actions.append(args))

    def fail() -> bool:
        raise RuntimeError("test network failure")

    monkeypatch.setattr(MODULE, "poll", fail)
    with pytest.raises(RuntimeError, match="network failure"):
        MODULE.main()
    assert actions[-1][0] == "outcome" and "failed" in actions[-1][1]


def test_main_exits_after_complete_return(monkeypatch: pytest.MonkeyPatch) -> None:
    """A completed cohort terminates the periodic monitor without another poll."""
    monkeypatch.setattr(MODULE, "poll", lambda: True)
    MODULE.main()
