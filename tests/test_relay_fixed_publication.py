"""Isolated publication tests for the fixed relay PMF watcher."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent.relay_umbrella_analysis_io import digest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/watch_relay_fixed_pmf.py"
SPEC = importlib.util.spec_from_file_location("watch_relay_fixed_pmf", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
WATCHER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WATCHER)


def setup_publication(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    checkpoint = tmp_path / "checkpoint_20ps"
    evidence = tmp_path / "evidence"
    plots = tmp_path / "plots"
    source = tmp_path / "source"
    for directory in (checkpoint, evidence, plots, source / "scripts"):
        directory.mkdir(parents=True, exist_ok=True)
    (checkpoint / "analysis.json").write_text(
        json.dumps({"checkpoint": {"actual_horizon_ps": 20.0}})
    )
    return checkpoint, evidence, plots, source


def write_report(checkpoint: Path, report: Path, source_sha: str | None = None) -> None:
    parent = checkpoint / "analysis.json"
    report.write_text(
        json.dumps(
            {"source_sha256": {str(parent.resolve()): source_sha or digest(parent)}}
        )
    )


def install_plot_stub(source: Path) -> None:
    plot_script = source / "scripts/plot_relay_fixed_pmf.py"
    plot_script.write_text(
        "from pathlib import Path\n"
        "def plot_report(report, output_dir):\n"
        "    output = Path(output_dir)\n"
        "    output.mkdir(parents=True, exist_ok=True)\n"
        "    pmf = output / 'checkpoint_pmf.png'\n"
        "    pmf.write_text('new pmf')\n"
        "    aux = output / 'checkpoint_overlap.png'\n"
        "    aux.write_text('overlap')\n"
        "    return {'plots': [str(pmf), str(aux)]}\n"
    )


def test_missing_retained_report_runs_fixed_analysis_and_publishes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkpoint, evidence, plots, source = setup_publication(tmp_path)
    install_plot_stub(source)
    calls: list[list[str]] = []

    def fake_run(command: list[str], *, check: bool) -> None:
        assert check is True
        calls.append(command)
        if "analyze_relay_fixed_pmf.py" in " ".join(command):
            report = Path(command[command.index("--output") + 1])
            write_report(checkpoint, report)

    copied_to: list[Path] = []
    original_copy2 = WATCHER.shutil.copy2

    def capture_copy2(src: str | Path, dst: str | Path) -> str:
        copied_to.append(Path(dst))
        return original_copy2(src, dst)

    monkeypatch.setattr(WATCHER.subprocess, "run", fake_run)
    monkeypatch.setattr(WATCHER.shutil, "copy2", capture_copy2)
    WATCHER.publish(checkpoint, evidence, plots, source, tmp_path)

    analyzer = next(
        call for call in calls if "analyze_relay_fixed_pmf.py" in " ".join(call)
    )
    assert analyzer[analyzer.index("--replicates") + 1] == "200"
    assert analyzer[analyzer.index("--workers") + 1] == "2"
    assert [call[call.index("--phase") + 1] for call in calls if "--phase" in call] == [
        "intent",
        "outcome",
    ]
    alias = plots / "pmf_latest.png"
    assert alias.read_text() == "new pmf"
    assert any(path.name == "pmf_latest.png.fixed_pending" for path in copied_to)
    assert not (plots / "pmf_latest.png.fixed_pending").exists()
    assert (plots / "checkpoint_overlap.png").exists()
    publication = json.loads((evidence / "fixed_20ps_publication.json").read_text())
    assert publication["parent_sha256"] == digest(checkpoint / "analysis.json")
    assert publication["native_jobs_changed"] is False


def test_retained_report_must_match_actual_parent_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkpoint, evidence, plots, source = setup_publication(tmp_path)
    report = evidence / "fixed_20ps_analysis.json"
    write_report(checkpoint, report, "0" * 64)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        WATCHER.subprocess, "run", lambda command, check: calls.append(command)
    )
    with pytest.raises(ValueError, match="another source report"):
        WATCHER.publish(checkpoint, evidence, plots, source, tmp_path)
    assert len(calls) == 1
    assert calls[0][calls[0].index("--phase") + 1] == "intent"


def test_changed_previous_publication_parent_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkpoint, evidence, plots, source = setup_publication(tmp_path)
    publication = evidence / "fixed_20ps_publication.json"
    publication.write_text(json.dumps({"parent_sha256": "wrong"}))
    called = False

    def fail_run(*args: Any, **kwargs: Any) -> None:
        nonlocal called
        called = True

    monkeypatch.setattr(WATCHER.subprocess, "run", fail_run)
    with pytest.raises(ValueError, match="Previously published parent changed"):
        WATCHER.publish(checkpoint, evidence, plots, source, tmp_path)
    assert not called


def test_valid_existing_publication_is_noop(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkpoint, evidence, plots, source = setup_publication(tmp_path)
    parent = checkpoint / "analysis.json"
    (evidence / "fixed_20ps_publication.json").write_text(
        json.dumps({"parent_sha256": digest(parent)})
    )
    monkeypatch.setattr(
        WATCHER.subprocess,
        "run",
        lambda *args, **kwargs: pytest.fail(
            "already published source should be a no-op"
        ),
    )
    WATCHER.publish(checkpoint, evidence, plots, source, tmp_path)


def test_main_rejects_invalid_runtime_before_creating_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence = tmp_path / "evidence"
    monkeypatch.setattr(
        "sys.argv",
        [
            str(SCRIPT),
            "--checkpoints",
            str(tmp_path),
            "--evidence",
            str(evidence),
            "--plots",
            str(tmp_path),
            "--source",
            str(tmp_path),
            "--journal-root",
            str(tmp_path),
            "--minimum-horizon",
            "20",
            "--max-hours",
            "0",
        ],
    )
    with pytest.raises(ValueError, match="Positive runtime"):
        WATCHER.main()
    assert not evidence.exists()


def test_main_uses_controlled_deadline_to_scan_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkpoints = tmp_path / "checkpoints"
    checkpoint = checkpoints / "checkpoint_20ps"
    checkpoint.mkdir(parents=True)
    (checkpoint / "analysis.json").write_text(
        json.dumps({"checkpoint": {"actual_horizon_ps": 20.0}})
    )
    evidence = tmp_path / "evidence"
    args = [
        str(SCRIPT),
        "--checkpoints",
        str(checkpoints),
        "--evidence",
        str(evidence),
        "--plots",
        str(tmp_path),
        "--source",
        str(tmp_path),
        "--journal-root",
        str(tmp_path),
        "--minimum-horizon",
        "20",
        "--max-hours",
        "0.001",
    ]
    monotonic_values = iter([0.0, 0.0, 5.0])
    published: list[Path] = []
    monkeypatch.setattr("sys.argv", args)
    monkeypatch.setattr(WATCHER.time, "monotonic", lambda: next(monotonic_values))
    monkeypatch.setattr(WATCHER.time, "sleep", lambda _: None)
    monkeypatch.setattr(WATCHER, "publish", lambda path, *args: published.append(path))
    WATCHER.main()
    assert published == [checkpoint]
    assert evidence.is_dir()
