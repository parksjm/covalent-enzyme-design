"""Focused tests for matched cumulative MBAR coverage selection."""

import json
from argparse import Namespace
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent.relay_mbar_progression import plot_progression, select_horizon
from scripts import watch_relay_umbrella_analysis as watcher


def test_select_horizon_requires_all_windows_and_caps_ahead_blocks(
    tmp_path: Path,
) -> None:
    windows = [{"label": "w1"}, {"label": "w2"}]
    receipts = []
    for label, durations in (("w1", [5, 5]), ("w2", [5, 5])):
        for index, duration in enumerate(durations, 1):
            stage = tmp_path / label / f"{index + 2:02d}_production_5ps"
            stage.mkdir(parents=True)
            path = stage / "analysis_acceptance.json"
            path.write_text(
                json.dumps(
                    {
                        "window": label,
                        "accepted": True,
                        "sampled_duration_ps": duration,
                    }
                )
            )
            receipts.append(path)

    selected = select_horizon(sorted(receipts), windows, 5)

    assert selected is not None
    paths, durations = selected
    assert len(paths) == 2
    assert durations == {"w1": 5, "w2": 5}


def test_select_horizon_waits_for_every_window(tmp_path: Path) -> None:
    stage = tmp_path / "w1" / "03_production_5ps"
    stage.mkdir(parents=True)
    path = stage / "analysis_acceptance.json"
    path.write_text(
        json.dumps(
            {
                "window": "w1",
                "accepted": True,
                "sampled_duration_ps": 5,
            }
        )
    )

    assert select_horizon([path], [{"label": "w1"}, {"label": "w2"}], 5) is None


def test_select_horizon_rejects_invalid_receipt(tmp_path: Path) -> None:
    stage = tmp_path / "w1" / "03_production_5ps"
    stage.mkdir(parents=True)
    path = stage / "analysis_acceptance.json"
    path.write_text(
        json.dumps(
            {
                "window": "w1",
                "accepted": False,
                "sampled_duration_ps": 5,
            }
        )
    )

    with pytest.raises(ValueError, match="Invalid accepted"):
        select_horizon([path], [{"label": "w1"}], 5)


def test_select_horizon_rejects_hidden_production_block(tmp_path: Path) -> None:
    paths = []
    for index in (3, 5):
        stage = tmp_path / "w1" / f"{index:02d}_production_5ps"
        stage.mkdir(parents=True)
        path = stage / "analysis_acceptance.json"
        path.write_text(
            json.dumps(
                {
                    "window": "w1",
                    "accepted": True,
                    "sampled_duration_ps": 5,
                }
            )
        )
        paths.append(path)

    with pytest.raises(ValueError, match="Noncontiguous"):
        select_horizon(paths, [{"label": "w1"}], 5)


def _connected_report(reference_value: float) -> dict[str, Any]:
    return {
        "status": "preliminary_connected_not_converged",
        "physical_hamiltonian_sha256": "same-hamiltonian",
        "temperature_K": 300,
        "path_3cv_sha256": "same-path",
        "projection": "same-projection",
        "es_search_bounds": [-0.1, 0.1],
        "acyl_search_bounds": [0.9, 1.1],
        "pmf": {
            "edges": [0.0, 0.5, 1.0],
            "es_search_bounds": [-0.1, 0.1],
            "acyl_search_bounds": [0.9, 1.1],
            "centers": [0.25, 0.75],
            "pmf_kcal_mol": [reference_value, reference_value + 2.0],
            "reference_s": 0.25,
        },
    }


def test_plot_progression_references_all_curves_to_same_es_bin(tmp_path: Path) -> None:
    output = tmp_path / "progression.png"

    result = plot_progression(
        [(5.0, _connected_report(0.0)), (10.0, _connected_report(1.0))], output
    )

    assert output.is_file()
    assert output.stat().st_size > 0
    assert result["common_reference_s"] == 0.25
    assert result["included_horizons_ps"] == [5.0, 10.0]
    assert result["reference_offsets_kcal_mol"] == [(5.0, 0.0), (10.0, 1.0)]
    assert result["uncertainty"] == "not calculated"


def test_plot_progression_marks_disconnected_only_progress_as_unavailable(
    tmp_path: Path,
) -> None:
    output = tmp_path / "disconnected.png"

    plot_progression(
        [(5.0, {"status": "disconnected_no_global_PMF", "pmf": None})], output
    )

    assert output.is_file()


def test_plot_progression_rejects_existing_output(tmp_path: Path) -> None:
    output = tmp_path / "existing.png"
    output.write_bytes(b"preserved")

    with pytest.raises(FileExistsError):
        plot_progression([], output)


def test_plot_progression_reports_curve_without_common_reference_support(
    tmp_path: Path,
) -> None:
    report = _connected_report(1.0)
    report["pmf"]["pmf_kcal_mol"][0] = None
    result = plot_progression(
        [(5.0, _connected_report(0.0)), (10.0, report)],
        tmp_path / "progression.png",
    )
    assert result["included_horizons_ps"] == [5.0]
    assert result["omitted_horizons_ps"] == [10.0]


def test_plot_progression_rejects_different_projection_metadata(
    tmp_path: Path,
) -> None:
    second = _connected_report(1.0)
    second["projection"] = "different-projection"
    with pytest.raises(ValueError, match="physics, projection"):
        plot_progression(
            [(5.0, _connected_report(0.0)), (10.0, second)],
            tmp_path / "mismatch.png",
        )


def test_watcher_emits_dynamic_horizons_once_and_resumes_previous(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    windows = [{"label": "w1"}, {"label": "w2"}]
    receipts_root = tmp_path / "receipts"
    receipts = []
    for label in ("w1", "w2"):
        for index in range(3, 10):
            stage = receipts_root / label / f"{index:02d}_production_5ps"
            stage.mkdir(parents=True)
            path = stage / "analysis_acceptance.json"
            path.write_text(
                json.dumps(
                    {
                        "window": label,
                        "accepted": True,
                        "sampled_duration_ps": 5,
                    }
                )
            )
            receipts.append(path)
    request = tmp_path / "request.json"
    request.write_text(
        json.dumps(
            {
                "windows": windows,
                "physical_hamiltonian_sha256": "same-hamiltonian",
                "es_search_bounds": [-0.1, 0.1],
            }
        )
    )
    output = tmp_path / "analysis"
    args = Namespace(
        request=request,
        receipts_root=receipts_root,
        receipt_pattern="analysis_acceptance.json",
        output_root=output,
        max_minutes=10080,
        poll_seconds=60,
        journal_root=tmp_path,
        resume_output=True,
        once=True,
    )
    calls = []
    comparisons = []

    def analyze(_request):
        calls.append(_request)
        return {
            "status": "preliminary_connected_not_converged",
            "physical_hamiltonian_sha256": "same-hamiltonian",
            "gates": {},
            "pmf": {
                "centers": [0.25],
                "pmf_kcal_mol": [0.0],
                "reference_s": 0.25,
            },
        }

    monkeypatch.setattr(watcher, "analyze_request", analyze)
    monkeypatch.setattr(
        watcher,
        "compare_reports",
        lambda *_: comparisons.append(True) or {"passed": None},
    )
    monkeypatch.setattr(
        watcher, "plot_report", lambda _, path: path.write_bytes(b"png")
    )
    monkeypatch.setattr(
        watcher, "plot_progression", lambda _, path: path.write_bytes(b"png")
    )
    monkeypatch.setattr(watcher, "journal", lambda *_args, **_kwargs: None)

    watcher.watch(args)

    assert len(calls) == 7
    assert len(comparisons) == 6
    assert (
        json.loads((output / "latest_checkpoint.json").read_text())["actual_horizon_ps"]
        == 35
    )
    assert (output / "checkpoint_35ps" / "source.json").is_file()

    calls.clear()
    comparisons.clear()
    watcher.watch(args)

    assert calls == []
    assert (
        json.loads((output / "latest_checkpoint.json").read_text())["actual_horizon_ps"]
        == 35
    )
