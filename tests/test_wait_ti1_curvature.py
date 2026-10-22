"""Safety tests for the conditional TI1 curvature controller."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from scripts import wait_nyl12_ti1_curvature as controller


def _write_source_snapshot(snapshot: Path) -> Path:
    """Create a one-file source snapshot with its recorded digest."""
    snapshot.mkdir()
    source = snapshot / "source.py"
    source.write_text("source is unchanged\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    (snapshot / "source_identity.json").write_text(json.dumps({"source.py": digest}))
    return source


def _write_completed_frequency(case: Path, exit_status: int = 0) -> None:
    """Write the required completed status record."""
    status = case / "frequency_controller_01" / "job" / "status.json"
    status.parent.mkdir(parents=True, exist_ok=True)
    status.write_text(
        json.dumps({"ended_utc": "2026-10-22T12:00:00Z", "exit_status": exit_status})
    )


@pytest.mark.parametrize("wait_hours,poll_seconds", [(0.0, 1.0), (1.0, 0.0)])
def test_wait_rejects_nonpositive_timing(
    tmp_path: Path, wait_hours: float, poll_seconds: float
) -> None:
    """Polling duration and deadline must both be positive."""
    with pytest.raises(ValueError):
        controller.wait_for_completion(
            tmp_path / "status.json", wait_hours, poll_seconds
        )


def test_verify_sources_rejects_empty_manifest(tmp_path: Path) -> None:
    """An empty identity manifest cannot establish a verified source set."""
    snapshot = tmp_path / "snapshot"
    snapshot.mkdir()
    (snapshot / "source_identity.json").write_text("{}")

    with pytest.raises(ValueError):
        controller.verify_sources(snapshot)


@pytest.mark.parametrize("relative_path", ["../outside.py", "sub/../../outside.py"])
def test_verify_sources_rejects_path_escape(tmp_path: Path, relative_path: str) -> None:
    """Identity paths must remain inside the source snapshot."""
    snapshot = tmp_path / "snapshot"
    snapshot.mkdir()
    (snapshot / "source_identity.json").write_text(
        json.dumps({relative_path: "0" * 64})
    )

    with pytest.raises(ValueError):
        controller.verify_sources(snapshot)


def test_verify_sources_rejects_changed_digest(tmp_path: Path) -> None:
    """A changed source is detected before any downstream work starts."""
    snapshot = tmp_path / "snapshot"
    source = _write_source_snapshot(snapshot)
    source.write_text("changed source\n")

    with pytest.raises(ValueError):
        controller.verify_sources(snapshot)


def test_wait_accepts_successfully_completed_job(tmp_path: Path) -> None:
    """A completed zero-exit frequency job satisfies the prerequisite."""
    status = tmp_path / "status.json"
    status.write_text(
        json.dumps({"ended_utc": "2026-10-22T12:00:00Z", "exit_status": 0})
    )

    controller.wait_for_completion(status, wait_hours=1.0, poll_seconds=1.0)


def test_wait_rejects_failed_completed_job(tmp_path: Path) -> None:
    """A completed nonzero frequency job blocks curvature analysis."""
    status = tmp_path / "status.json"
    status.write_text(
        json.dumps({"ended_utc": "2026-10-22T12:00:00Z", "exit_status": 7})
    )

    with pytest.raises(RuntimeError):
        controller.wait_for_completion(status, wait_hours=1.0, poll_seconds=1.0)


def test_wait_times_out_if_status_never_appears(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An absent status record eventually stops the controller."""
    clock_values = iter([0.0, 2.0])
    monkeypatch.setattr(controller.time, "monotonic", lambda: next(clock_values))
    monkeypatch.setattr(controller.time, "sleep", lambda seconds: None)

    with pytest.raises(TimeoutError):
        controller.wait_for_completion(
            tmp_path / "missing.json", wait_hours=1e-6, poll_seconds=1.0
        )


def _install_run_prerequisites(
    tmp_path: Path,
) -> tuple[Path, Path, Path, Path]:
    """Create valid snapshot, case status, and journal paths."""
    snapshot = tmp_path / "snapshot"
    _write_source_snapshot(snapshot)
    case = tmp_path / "case"
    case.mkdir()
    _write_completed_frequency(case)
    journal = tmp_path / "actions.jsonl"
    return case, snapshot, journal, snapshot / "source.py"


def test_run_analysis_calls_analyzer_with_explicit_settings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A verified completed job launches analysis with recorded settings."""
    case, snapshot, journal, _source = _install_run_prerequisites(tmp_path)
    calls: list[tuple[tuple[Any, ...], dict[str, Any]]] = []

    def fake_run(*args: Any, **kwargs: Any) -> SimpleNamespace:
        calls.append((args, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(controller.subprocess, "run", fake_run)

    result = controller.run_analysis(
        case, snapshot, journal, wait_hours=1.0, poll_seconds=1.0
    )

    assert result == 0
    assert len(calls) == 3  # intent record, analyzer, outcome record
    analyzer_call = next(
        args
        for args, _kwargs in calls
        if any("analyze_nyl12_ti1_curvature" in str(part) for part in args[0])
    )
    command = list(analyzer_call[0])
    assert "--frequency-dir" in command
    assert str(case / "frequency_01") == command[command.index("--frequency-dir") + 1]
    assert "--output" in command
    assert (
        str(case / "conditional_curvature_01") == command[command.index("--output") + 1]
    )
    expected = {
        "--step-bohr": "0.0002",
        "--stationarity-tolerance": "0.00003",
        "--curvature-tolerance": "0.000001",
        "--derivative-relative-tolerance": "0.001",
    }
    for option, value in expected.items():
        assert command[command.index(option) + 1] == value
    assert calls[0][1]["check"] is True
    assert calls[1][1]["check"] is False
    assert calls[2][1]["check"] is True


@pytest.mark.parametrize("failure", ["frequency", "hash"])
def test_failed_prerequisite_or_hash_never_invokes_analyzer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
) -> None:
    """A failed frequency result or changed source blocks analyzer execution."""
    case, snapshot, journal, source = _install_run_prerequisites(tmp_path)
    if failure == "frequency":
        _write_completed_frequency(case, exit_status=3)
    else:
        source.write_text("mutated source\n")
    calls: list[tuple[Any, ...]] = []

    def fake_run(*args: Any, **kwargs: Any) -> SimpleNamespace:
        calls.append(args)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(controller.subprocess, "run", fake_run)

    with pytest.raises((RuntimeError, ValueError)):
        controller.run_analysis(
            case, snapshot, journal, wait_hours=1.0, poll_seconds=1.0
        )

    assert not any(
        any("analyze_nyl12_ti1_curvature" in str(part) for part in args[0])
        for args in calls
    )


def test_run_analysis_refuses_existing_output_before_analyzer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Existing assessment output is preserved and prevents analyzer launch."""
    case, snapshot, journal, _source = _install_run_prerequisites(tmp_path)
    output = case / "conditional_curvature_01"
    output.mkdir()
    calls: list[tuple[Any, ...]] = []
    monkeypatch.setattr(
        controller.subprocess,
        "run",
        lambda *args, **kwargs: calls.append(args) or SimpleNamespace(returncode=0),
    )

    with pytest.raises((FileExistsError, ValueError)):
        controller.run_analysis(
            case, snapshot, journal, wait_hours=1.0, poll_seconds=1.0
        )

    assert not any(
        any("analyze_nyl12_ti1_curvature" in str(part) for part in args[0])
        for args in calls
    )


def test_main_requires_explicit_controller_settings(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """CLI requires all case, snapshot, journal, and polling inputs."""
    monkeypatch.setattr("sys.argv", ["wait_nyl12_ti1_curvature"])

    with pytest.raises(SystemExit) as error:
        controller.main()

    assert error.value.code == 2
