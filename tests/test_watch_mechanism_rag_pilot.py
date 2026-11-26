"""Read-only lifecycle tests for the mechanism RAG monitor."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Any

import pytest

from scripts import watch_mechanism_rag_pilot as monitor


def _write_status(path: Path, *, ended: str | None, exit_status: int | None) -> None:
    path.write_text(json.dumps({"ended_utc": ended, "exit_status": exit_status}))


def _summary(
    *,
    response_count: int = 20,
    citation_valid_count: int = 20,
    issues: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "expected_calls": 20,
        "response_count": response_count,
        "citation_valid_count": citation_valid_count,
        "issues": [] if issues is None else issues,
    }


def test_terminal_success_writes_report_without_duplicate_notification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    controller = tmp_path / "controller.json"
    _write_status(controller, ended="done", exit_status=0)
    notices: list[tuple[str, Path, str]] = []
    records: list[str] = []
    monkeypatch.setattr(monitor, "notify", lambda *args: notices.append(args))
    monkeypatch.setattr(monitor, "record", lambda phase, detail: records.append(phase))
    summary = _summary()

    monitor.monitor(
        tmp_path / "pilot",
        controller,
        tmp_path / "out",
        lambda *_: summary,
        lambda value: f"{value['response_count']}/{value['expected_calls']} responses",
        monotonic=lambda: 0.0,
    )

    report = json.loads((tmp_path / "out/final_report.json").read_text())
    assert report["summary"] == summary
    assert report["monitor_deadline_reached"] is False
    assert (tmp_path / "out/MORNING_REPORT.md").is_file()
    assert notices == []
    assert records.count("outcome") == 2


def test_incomplete_error_sends_one_truthful_followup_with_actual_counts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    controller = tmp_path / "controller.json"
    _write_status(controller, ended="failed", exit_status=1)
    summary = _summary(
        response_count=13, citation_valid_count=11, issues=["missing response"]
    )
    notices: list[tuple[str, Path, str]] = []
    monkeypatch.setattr(monitor, "record", lambda *_: None)
    monkeypatch.setattr(monitor, "notify", lambda *args: notices.append(args))

    monitor.monitor(
        tmp_path / "pilot",
        controller,
        tmp_path / "out",
        lambda *_: summary,
        lambda value: (
            f"responses={value['response_count']}/{value['expected_calls']}; "
            f"citation-valid={value['citation_valid_count']}/{value['expected_calls']}"
        ),
        monotonic=lambda: 0.0,
    )

    assert len(notices) == 1
    assert notices[0][0] == "error"
    assert "responses=13/20" in notices[0][2]
    assert "citation-valid=11/20" in notices[0][2]


def test_running_controller_is_polled_until_terminal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    controller = tmp_path / "controller.json"
    _write_status(controller, ended=None, exit_status=None)
    polls = 0
    sleeps: list[float] = []
    now = 0.0
    notices: list[tuple[Any, ...]] = []
    monkeypatch.setattr(monitor, "record", lambda *_: None)
    monkeypatch.setattr(monitor, "notify", lambda *args: notices.append(args))

    def summarize(*_: Path) -> dict[str, Any]:
        nonlocal polls
        polls += 1
        return _summary()

    def sleep(seconds: float) -> None:
        nonlocal now
        sleeps.append(seconds)
        now += seconds
        _write_status(controller, ended="done", exit_status=0)

    monitor.monitor(
        tmp_path / "pilot",
        controller,
        tmp_path / "out",
        summarize,
        lambda _: "complete",
        hours=1,
        interval=7,
        sleep=sleep,
        monotonic=lambda: now,
    )

    assert polls == 2
    assert sleeps == [7]
    assert notices == []
    assert (
        json.loads((tmp_path / "out/final_report.json").read_text())[
            "native_exit_status"
        ]
        == 0
    )


def test_deadline_sends_question_and_leaves_controller_unchanged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    controller = tmp_path / "controller.json"
    _write_status(controller, ended=None, exit_status=None)
    original = controller.read_bytes()
    ticks = iter((0.0, 4.0))
    notices: list[tuple[str, Path, str]] = []
    records: list[tuple[str, str]] = []
    monkeypatch.setattr(monitor, "record", lambda *args: records.append(args))
    monkeypatch.setattr(monitor, "notify", lambda *args: notices.append(args))

    monitor.monitor(
        tmp_path / "pilot",
        controller,
        tmp_path / "out",
        lambda *_: _summary(response_count=3, citation_valid_count=2),
        lambda value: f"actual responses={value['response_count']}",
        hours=1 / 3600,
        interval=2,
        sleep=lambda _: pytest.fail("deadline path must not sleep"),
        monotonic=lambda: next(ticks),
    )

    assert controller.read_bytes() == original
    assert len(notices) == 1
    assert notices[0][0] == "question"
    assert "Audit deadline reached" in notices[0][2]
    assert "actual responses=3" in notices[0][2]
    assert (
        json.loads((tmp_path / "out/final_report.json").read_text())[
            "monitor_deadline_reached"
        ]
        is True
    )
    assert any(phase == "intent" and "follow-up" in detail for phase, detail in records)


@pytest.mark.parametrize(("hours", "interval"), [(0, 1), (1, 0.5), (-1, 60)])
def test_invalid_bounds_fail_before_creating_output_or_reading_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    hours: float,
    interval: float,
) -> None:
    output = tmp_path / "out"
    monkeypatch.setattr(monitor, "record", lambda *_: pytest.fail("must not journal"))

    with pytest.raises(ValueError, match="Positive duration"):
        monitor.monitor(
            tmp_path / "pilot",
            tmp_path / "missing-controller.json",
            output,
            lambda *_: pytest.fail("must not summarize"),
            lambda _: "unused",
            hours=hours,
            interval=interval,
        )

    assert not output.exists()


def test_notify_uses_the_trusted_existing_config_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    agent_root = tmp_path / "agent"
    calls: list[tuple[list[str], dict[str, Any]]] = []

    def fake_run(command: list[str], **kwargs: Any) -> SimpleNamespace:
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(monitor, "AGENT_ROOT", agent_root)
    monkeypatch.setattr(monitor.subprocess, "run", fake_run)

    monitor.notify("question", tmp_path / "audit", "actual counts only")

    command, kwargs = calls[0]
    assert command[0] == str(agent_root / ".venv/bin/mechanism-agent")
    assert command[command.index("--notifications-config") + 1] == str(
        agent_root / "deployment/notifications.toml"
    )
    assert command[command.index("--event") + 1] == "question"
    assert kwargs["check"] is True


def test_record_uses_both_journals_and_stops_after_journal_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(monitor, "SCIENCE_ROOT", tmp_path / "science")
    monkeypatch.setattr(monitor, "AGENT_ROOT", tmp_path / "agent")
    calls: list[tuple[list[str], dict[str, Any]]] = []

    def fail_first(command: list[str], **kwargs: Any) -> SimpleNamespace:
        calls.append((command, kwargs))
        raise OSError("journal unavailable")

    monkeypatch.setattr(monitor.subprocess, "run", fail_first)
    with pytest.raises(OSError, match="journal unavailable"):
        monitor.record("intent", "test")

    assert len(calls) == 1
    assert calls[0][0][calls[0][0].index("--log") + 1] == str(
        tmp_path / "science/output_agent_actions/20261021.jsonl"
    )
    assert calls[0][1]["check"] is True


def test_journal_failure_prevents_summary_and_notification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    controller = tmp_path / "controller.json"
    _write_status(controller, ended="done", exit_status=1)
    notices: list[tuple[Any, ...]] = []

    def fail_on_notify_intent(phase: str, detail: str) -> None:
        if "Send one factual follow-up" in detail:
            raise OSError("journal unavailable")

    monkeypatch.setattr(monitor, "record", fail_on_notify_intent)
    monkeypatch.setattr(monitor, "notify", lambda *args: notices.append(args))

    with pytest.raises(OSError, match="journal unavailable"):
        monitor.monitor(
            tmp_path / "pilot",
            controller,
            tmp_path / "out",
            lambda *_: _summary(response_count=1),
            lambda _: "follow-up",
            monotonic=lambda: 0.0,
        )

    assert notices == []
    assert (tmp_path / "out/final_report.json").is_file()


def test_main_uses_frozen_sibling_reporting_helper(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: list[tuple[Any, ...]] = []
    reporting = SimpleNamespace(
        summarize_pilot=lambda *_: {},
        format_notification=lambda _: "summary",
    )
    monkeypatch.setitem(sys.modules, "reporting", reporting)
    monkeypatch.setattr(
        monitor, "monitor", lambda *args, **kwargs: captured.append((args, kwargs))
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "watch_mechanism_rag_pilot.py",
            "--pilot",
            str(tmp_path / "pilot"),
            "--controller-status",
            str(tmp_path / "controller.json"),
            "--output",
            str(tmp_path / "out"),
            "--hours",
            "0.25",
            "--interval",
            "5",
        ],
    )

    monitor.main()

    args, kwargs = captured[0]
    assert args[:3] == (
        tmp_path / "pilot",
        tmp_path / "controller.json",
        tmp_path / "out",
    )
    assert args[3] is reporting.summarize_pilot
    assert args[4] is reporting.format_notification
    assert kwargs == {"hours": 0.25, "interval": 5.0}


def test_main_reports_monitor_failure_truthfully(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    notices: list[tuple[str, Path, str]] = []
    records: list[tuple[str, str]] = []
    monkeypatch.setitem(
        sys.modules,
        "reporting",
        SimpleNamespace(
            summarize_pilot=lambda *_: {}, format_notification=lambda _: "unused"
        ),
    )
    monkeypatch.setattr(
        monitor,
        "monitor",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("report failed")),
    )
    monkeypatch.setattr(monitor, "record", lambda *args: records.append(args))
    monkeypatch.setattr(monitor, "notify", lambda *args: notices.append(args))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "watch_mechanism_rag_pilot.py",
            "--pilot",
            str(tmp_path / "pilot"),
            "--controller-status",
            str(tmp_path / "controller.json"),
            "--output",
            str(tmp_path / "out"),
        ],
    )

    with pytest.raises(RuntimeError, match="report failed"):
        monitor.main()

    assert len(notices) == 1
    assert notices[0][0] == "error"
    assert "No claim about completed answer count" in notices[0][2]
    assert "20 calls completed" not in notices[0][2]
    assert records == [
        ("intent", "Report an audit failure without changing the pilot"),
        ("outcome", "Audit failure notification command completed"),
    ]


def test_main_does_not_notify_if_error_journal_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    notices: list[tuple[Any, ...]] = []
    monkeypatch.setitem(
        sys.modules,
        "reporting",
        SimpleNamespace(
            summarize_pilot=lambda *_: {}, format_notification=lambda _: "unused"
        ),
    )
    monkeypatch.setattr(
        monitor,
        "monitor",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("summary failed")),
    )
    monkeypatch.setattr(
        monitor, "record", lambda *_: (_ for _ in ()).throw(OSError("journal down"))
    )
    monkeypatch.setattr(monitor, "notify", lambda *args: notices.append(args))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "watch_mechanism_rag_pilot.py",
            "--pilot",
            str(tmp_path / "pilot"),
            "--controller-status",
            str(tmp_path / "controller.json"),
            "--output",
            str(tmp_path / "out"),
        ],
    )

    with pytest.raises(OSError, match="journal down"):
        monitor.main()

    assert notices == []
