"""Focused tests for safe watcher evidence reuse and latest aliases."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest


SPEC = importlib.util.spec_from_file_location(
    "watch_direct_fixed_pmf",
    Path(__file__).parents[1] / "scripts/watch_direct_fixed_pmf.py",
)
assert SPEC is not None and SPEC.loader is not None
WATCHER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WATCHER)
from nyl12_covalent import direct_fixed_watch_helpers as HELPERS  # noqa: E402


def _sha(path: Path) -> str:
    return WATCHER.digest(path)


def test_saved_report_reuse_verifies_all_source_hashes(tmp_path: Path) -> None:
    source = tmp_path / "request.json"
    source.write_text("{}")
    report = tmp_path / "direct_35ps_analysis.json"
    report.write_text(
        json.dumps(
            {
                "horizon_ps": 35,
                "windows": 63,
                "protocol": "fixed5ps_once_per_continuous_production_window_all_frames",
                "source_sha256": {str(source): _sha(source)},
            }
        )
    )
    assert WATCHER._verify_report(report, 35)["windows"] == 63
    source.write_text("changed")
    with pytest.raises(ValueError, match="source changed"):
        WATCHER._verify_report(report, 35)


def test_latest_alias_update_preserves_old_target(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    plots.mkdir()
    old_png, old_pdf = plots / "old.png", plots / "old.pdf"
    old_png.write_bytes(b"prior png")
    old_pdf.write_bytes(b"prior pdf")
    (plots / "pmf_latest.png").symlink_to(old_png.name)
    (plots / "pmf_latest.pdf").symlink_to(old_pdf.name)
    new_png, new_pdf = plots / "new.png", plots / "new.pdf"
    new_png.write_bytes(b"new png")
    new_pdf.write_bytes(b"new pdf")

    aliases = WATCHER._replace_latest(plots, {"plots": [str(new_png), str(new_pdf)]})

    assert (plots / "pmf_latest.png").resolve() == new_png
    assert (plots / "pmf_latest.pdf").resolve() == new_pdf
    assert aliases[str(plots / "pmf_latest.png")]["previous_sha256"] == _sha(old_png)
    assert old_png.read_bytes() == b"prior png"


def test_latest_alias_refuses_regular_file(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    plots.mkdir()
    (plots / "pmf_latest.png").write_bytes(b"do not replace")
    png, pdf = plots / "new.png", plots / "new.pdf"
    png.write_bytes(b"png")
    pdf.write_bytes(b"pdf")

    with pytest.raises(ValueError, match="non-symlink"):
        WATCHER._replace_latest(plots, {"plots": [str(png), str(pdf)]})


def test_matching_relay_report_uses_57_window_contract(tmp_path: Path) -> None:
    """Matching sampling horizons need not have identical window counts."""
    source = tmp_path / "request.json"
    source.write_text("{}")
    report = tmp_path / "fixed_35ps_analysis.json"
    report.write_text(
        json.dumps(
            {
                "horizon_ps": 35,
                "windows": 57,
                "protocol": "fixed5ps_once_per_continuous_production_window_all_frames",
                "source_sha256": {str(source): _sha(source)},
            }
        )
    )
    assert WATCHER._verify_report(report, 35, expected_windows=57)["windows"] == 57
    with pytest.raises(ValueError, match="wrong fixed-analysis identity"):
        WATCHER._verify_report(report, 35, expected_windows=63)


def _fixed_report(path: Path, horizon: int, windows: int, source: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "horizon_ps": horizon,
                "windows": windows,
                "protocol": "fixed5ps_once_per_continuous_production_window_all_frames",
                "source_sha256": {str(source): _sha(source)},
            }
        )
    )


def test_watch_waits_for_relay_then_reuses_direct_report_and_publishes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = tmp_path / "campaign"
    analysis = campaign / "analysis/fixed_burnin_20261126"
    analysis.mkdir(parents=True)
    request = analysis / "analysis_request.json"
    request.write_text("{}")
    _fixed_report(analysis / "direct_30ps_analysis.json", 30, 63, request)
    plan = tmp_path / "plan.json"
    plan.write_text("{}")
    plan_digest = _sha(plan)
    augmented = tmp_path / "augmented"
    augmented.mkdir()
    evidence = tmp_path / "evidence"
    relay_evidence = tmp_path / "relay_evidence"
    relay_evidence.mkdir()
    plots = tmp_path / "plots"
    source = tmp_path / "analysis_source"
    (source / "scripts").mkdir(parents=True)
    (source / "scripts/analyze_direct_fixed_pmf.py").write_text("")
    (source / "scripts/plot_direct_fixed_pmf.py").write_text("")
    relay_source = tmp_path / "relay_source.json"
    relay_source.write_text("relay source")
    direct_path = evidence / "direct_60ps_analysis.json"
    relay_path = relay_evidence / "fixed_60ps_analysis.json"
    monkeypatch.setattr(WATCHER, "HORIZONS_PS", (60.0,))
    monkeypatch.setattr(WATCHER, "_repository_root", lambda _campaign: tmp_path)
    monkeypatch.setattr(
        WATCHER,
        "_continuation_ready",
        lambda *_args: (True, {}, {str(plan): plan_digest}, augmented),
    )
    monkeypatch.setattr(WATCHER, "_record_action", lambda *_args: None)
    analyzed: list[Path] = []
    plotted: list[Path] = []

    def fake_child(
        _root: Path,
        _evidence: Path,
        action: str,
        _reason: str,
        command: list[str],
        _log: Path,
        _output: Path,
        _timeout: float,
        _source: Path,
    ) -> None:
        if action.endswith("analysis"):
            output = Path(command[command.index("--output") + 1])
            analyzed.append(output)
            _fixed_report(output, 60, 63, plan)
            return
        output = Path(command[command.index("--output-dir") + 1])
        plotted.append(output)
        png, pdf = output / "comparison.png", output / "comparison.pdf"
        png.write_bytes(b"png bytes")
        pdf.write_bytes(b"pdf bytes")
        metadata = {
            "direct_report_sha256": _sha(direct_path),
            "relay_report_sha256": _sha(relay_path),
            "horizons_matched": True,
            "plots": [str(png), str(pdf)],
        }
        (output / "comparison_60ps_vs_relay_60ps_test_source.json").write_text(
            json.dumps(metadata)
        )

    monkeypatch.setattr(WATCHER, "_journaled_child", fake_child)

    class StopAfterWait(Exception):
        pass

    monkeypatch.setattr(
        WATCHER.time, "sleep", lambda _seconds: (_ for _ in ()).throw(StopAfterWait())
    )
    args = SimpleNamespace(
        campaign=campaign,
        continuation_plan=plan,
        relay_evidence=relay_evidence,
        source=source,
        evidence=evidence,
        plots=plots,
        max_hours=0.01,
        replicates=200,
        workers=6,
    )
    with pytest.raises(StopAfterWait):
        WATCHER.watch(args)
    status = json.loads((evidence / "direct_fixed_pmf_watch_status.json").read_text())
    assert status["status"] == "waiting_for_matching_relay_report"
    assert direct_path.is_file() and len(analyzed) == 1 and not plotted

    _fixed_report(relay_path, 60, 57, relay_source)
    WATCHER.watch(args)

    assert analyzed == [direct_path]
    assert plotted == [plots]
    assert (evidence / "direct_60ps_publication.json").is_file()
    assert (plots / "pmf_latest.png").resolve() == plots / "comparison.png"
    assert (plots / "pmf_latest.pdf").resolve() == plots / "comparison.pdf"


def test_flat_plot_metadata_filters_horizon_and_checks_source_hashes(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.json"
    source.write_text("inputs")
    direct = tmp_path / "direct.json"
    relay = tmp_path / "relay.json"
    _fixed_report(direct, 60, 63, source)
    _fixed_report(relay, 60, 57, source)
    other = tmp_path / "comparison_55ps_vs_relay_55ps_old_source.json"
    other.write_text("{}")
    assert WATCHER._plot_metadata(tmp_path, direct, relay) is None

    png, pdf = tmp_path / "matched.png", tmp_path / "matched.pdf"
    png.write_bytes(b"png")
    pdf.write_bytes(b"pdf")
    metadata_path = tmp_path / "comparison_60ps_vs_relay_60ps_current_source.json"
    metadata_path.write_text(
        json.dumps(
            {
                "direct_report_sha256": _sha(direct),
                "relay_report_sha256": _sha(relay),
                "horizons_matched": True,
                "plots": [str(png), str(pdf)],
            }
        )
    )
    assert WATCHER._plot_metadata(tmp_path, direct, relay)["plots"] == [
        str(png),
        str(pdf),
    ]
    direct.write_text('{"horizon_ps":60,"changed":true}')
    with pytest.raises(ValueError, match="stale"):
        WATCHER._plot_metadata(tmp_path, direct, relay)


def test_journaled_child_sets_frozen_pythonpath_and_retains_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[str, str]] = []
    monkeypatch.setattr(
        HELPERS,
        "_record_action",
        lambda _root, phase, action, *_args: calls.append((phase, action)),
    )
    evidence = tmp_path / "evidence"
    source_root = tmp_path / "snapshot"
    (source_root / "src").mkdir(parents=True)
    marker = tmp_path / "child_pythonpath.txt"
    log = evidence / "watcher_logs/child.log"
    command = [
        __import__("sys").executable,
        "-c",
        "import os,sys; open(sys.argv[1], 'w').write(os.environ['PYTHONPATH']); sys.exit(7)",
        str(marker),
    ]
    with pytest.raises(subprocess.CalledProcessError) as error:
        HELPERS._journaled_child(
            tmp_path,
            evidence,
            "child test",
            "test failure handling",
            command,
            log,
            tmp_path / "unused-output.json",
            10.0,
            source_root,
        )

    assert error.value.returncode == 7
    assert str(source_root / "src") in marker.read_text().split(":")
    failure = json.loads((evidence / "child_test_failure.json").read_text())
    assert failure["returncode"] == 7 and log.is_file()
    assert calls == [("intent", "child test"), ("outcome", "child test")]
