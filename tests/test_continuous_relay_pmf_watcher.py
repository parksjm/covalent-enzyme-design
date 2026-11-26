"""Check checkpoint identity and derived-only publication by the new watcher."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SPEC = importlib.util.spec_from_file_location(
    "continuous_watcher",
    Path(__file__).resolve().parents[1] / "scripts/watch_continuous_relay_pmf.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_derived_publication_and_no_overwrite_of_changed_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    checkpoint, evidence, plots = [
        tmp_path / name for name in ("checkpoint_20ps", "evidence", "plots")
    ]
    checkpoint.mkdir()
    plots.mkdir()
    report = checkpoint / "analysis.json"
    report.write_text(json.dumps(dict(checkpoint=dict(actual_horizon_ps=20))))
    (checkpoint / "analysis_request.json").write_text("{}")
    calls = []

    def run(command: list[str], *, check: bool) -> None:
        assert check
        calls.append(command)
        script = Path(command[1]).name
        if script == "analyze_relay_continuous_pmf.py":
            prefix = Path(command[command.index("--output-stem") + 1])
            for suffix in ("analysis", "uncertainty"):
                prefix.with_name(prefix.name + f"_{suffix}.json").write_text("{}")
        elif script == "analyze_relay_pmf_uncertainty.py":
            Path(command[command.index("--output") + 1]).write_text("{}")
        elif script == "plot_relay_pmf_diagnostics.py":
            assert command[-2:] == ["--interval-method", "analytic"]
            stem = command[command.index("--stem") + 1]
            for suffix in ("pmf", "mbar_overlap", "neighbor_overlap"):
                for extension in ("png", "pdf"):
                    (plots / f"{stem}_{suffix}.{extension}").write_text("figure")

    monkeypatch.setattr(MODULE.subprocess, "run", run)
    MODULE.publish(checkpoint, evidence, plots, tmp_path, tmp_path)
    assert len(calls) == 6
    assert (plots / "pmf_latest.png").read_text() == "figure"
    MODULE.publish(checkpoint, evidence, plots, tmp_path, tmp_path)
    assert len(calls) == 6
    report.write_text(report.read_text() + " ")
    with pytest.raises(ValueError, match="checkpoint changed"):
        MODULE.publish(checkpoint, evidence, plots, tmp_path, tmp_path)


def test_partial_analysis_fails_before_plotting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkpoint = tmp_path / "checkpoint"
    checkpoint.mkdir()
    (checkpoint / "analysis.json").write_text(
        json.dumps(dict(checkpoint=dict(actual_horizon_ps=20)))
    )
    (checkpoint / "analysis_request.json").write_text("{}")
    (tmp_path / "continuous_20ps_analysis.json").write_text("{}")
    monkeypatch.setattr(MODULE.subprocess, "run", lambda *a, **kw: None)
    with pytest.raises(ValueError, match="Partial"):
        MODULE.publish(checkpoint, tmp_path, tmp_path, tmp_path, tmp_path)


def test_watcher_waits_for_complete_source_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    for name, horizon in [("05", 5), ("20", 20), ("25", 25)]:
        path = tmp_path / f"checkpoint_{name}ps"
        path.mkdir()
        (path / "analysis.json").write_text(
            json.dumps(dict(checkpoint=dict(actual_horizon_ps=horizon)))
        )
        if horizon != 25:
            (path / "analysis_request.json").write_text("{}")
    args = SimpleNamespace(
        checkpoints=tmp_path,
        evidence=tmp_path,
        plots=tmp_path,
        source=tmp_path,
        journal_root=tmp_path,
        minimum_horizon=20,
        max_hours=72,
    )
    monkeypatch.setattr(MODULE.argparse.ArgumentParser, "parse_args", lambda _: args)
    clock = iter([0, 1, 1000000])
    monkeypatch.setattr(MODULE.time, "monotonic", lambda: next(clock))
    monkeypatch.setattr(MODULE.time, "sleep", lambda _: None)
    seen = []
    monkeypatch.setattr(MODULE, "publish", lambda *a: seen.append(a[0].name))
    MODULE.main()
    assert seen == ["checkpoint_20ps"]
    args.max_hours = 0
    with pytest.raises(ValueError, match="Positive"):
        MODULE.main()
