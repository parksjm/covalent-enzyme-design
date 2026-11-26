"""Exercise derived checkpoint publication without launching scientific jobs."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SPEC = importlib.util.spec_from_file_location(
    "figure_watcher",
    Path(__file__).resolve().parents[1] / "scripts/watch_relay_pmf_figures.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_publishes_latest_figures_and_reuses_only_identical_checkpoint(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    checkpoint = tmp_path / "checkpoint_20ps"
    checkpoint.mkdir()
    report = checkpoint / "analysis.json"
    report.write_text(json.dumps(dict(checkpoint=dict(actual_horizon_ps=20))))
    (checkpoint / "analysis_request.json").write_text("{}")
    evidence, plots = tmp_path / "evidence", tmp_path / "plots"
    calls = []

    def run(command: list[str], *, check: bool) -> None:
        assert check
        calls.append(command)
        if any("analyze_relay_pmf_uncertainty.py" in str(item) for item in command):
            Path(command[command.index("--output") + 1]).write_text("{}")
        if any("plot_relay_pmf_diagnostics.py" in str(item) for item in command):
            assert command[-2:] == ["--interval-method", "block-bootstrap"]
            stem = command[command.index("--stem") + 1]
            plots.mkdir()
            for suffix in ("pmf", "mbar_overlap", "neighbor_overlap"):
                for extension in ("png", "pdf"):
                    (plots / f"{stem}_{suffix}.{extension}").write_text("figure")

    monkeypatch.setattr(MODULE.subprocess, "run", run)
    options = dict(replicates=200, workers=2)
    result = MODULE.publish_checkpoint(
        checkpoint, evidence, plots, tmp_path, tmp_path, **options
    )
    assert result.is_file()
    assert len(calls) == 4
    assert (plots / "pmf_latest.png").read_text() == "figure"
    assert (
        MODULE.publish_checkpoint(
            checkpoint, evidence, plots, tmp_path, tmp_path, **options
        )
        == result
    )
    assert len(calls) == 4
    report.write_text(report.read_text() + " ")
    with pytest.raises(ValueError, match="checkpoint changed"):
        MODULE.publish_checkpoint(
            checkpoint, evidence, plots, tmp_path, tmp_path, **options
        )


def test_bounded_watcher_selects_completed_horizons_only(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = []
    for label, horizon in [("05", 5), ("20", 20), ("25", 25)]:
        checkpoint = tmp_path / f"checkpoint_{label}ps"
        checkpoint.mkdir()
        (checkpoint / "analysis.json").write_text(
            json.dumps(dict(checkpoint=dict(actual_horizon_ps=horizon)))
        )
        if horizon != 25:
            (checkpoint / "analysis_request.json").write_text("{}")
    options = SimpleNamespace(
        checkpoints=tmp_path,
        evidence=tmp_path,
        plots=tmp_path,
        source_root=tmp_path,
        journal_root=tmp_path,
        minimum_horizon=20,
        replicates=200,
        workers=2,
        max_hours=1,
        poll_seconds=60,
    )
    monkeypatch.setattr(MODULE.argparse.ArgumentParser, "parse_args", lambda _: options)
    moments = iter([0, 0, 99999])
    monkeypatch.setattr(MODULE.time, "monotonic", lambda: next(moments))
    monkeypatch.setattr(MODULE.time, "sleep", lambda _: None)
    monkeypatch.setattr(
        MODULE, "publish_checkpoint", lambda *a, **kw: calls.append(a[0])
    )
    MODULE.main()
    assert calls == [tmp_path / "checkpoint_20ps"]
    options.poll_seconds = 1
    with pytest.raises(ValueError, match="runtime bound"):
        MODULE.main()
