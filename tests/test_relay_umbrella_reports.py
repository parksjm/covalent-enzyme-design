"""Report, bootstrap and readiness watcher tests."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent import relay_umbrella_bootstrap as bootstrap_module
from nyl12_covalent.relay_umbrella_analysis_reports import (
    chemical_diagnostics,
    compare_reports,
    plot_report,
)


def report() -> dict:
    """Provide a small supported PMF with explicitly resolved basin references."""
    return dict(
        physical_hamiltonian_sha256="a" * 64,
        temperature_K=300,
        projection="polyline",
        path_3cv_sha256="b" * 64,
        path_length_A=5.0,
        status="preliminary",
        N_k=[100, 120],
        overlap=[[0.8, 0.2], [0.2, 0.8]],
        pmf=dict(
            edges=[0, 0.25, 0.5, 0.75, 1.0],
            centers=[0.125, 0.375, 0.625, 0.875],
            es_search_bounds=[0, 0.25],
            acyl_search_bounds=[0.75, 1.0],
            pmf_kcal_mol=[0.0, 1.0, None, 0.5],
            minimum_to_barrier_kcal_mol=1.0,
            basin_minimum_difference_kcal_mol=0.5,
            reference_s=0.125,
            acyl_reference_s=0.875,
            barrier_s=0.375,
        ),
    )


def test_report_comparison_scope_and_unresolved() -> None:
    first, second = report(), report()
    assert compare_reports(second, first)["passed"]
    second["pmf"]["minimum_to_barrier_kcal_mol"] = 1.8
    assert not compare_reports(second, first)["passed"]
    second["pmf"]["minimum_to_barrier_kcal_mol"] = None
    assert compare_reports(second, first)["passed"] is None
    second["temperature_K"] = 301
    with pytest.raises(ValueError, match="temperature"):
        compare_reports(second, first)
    second = report()
    second["pmf"]["edges"] = [0, 1]
    with pytest.raises(ValueError, match="edges"):
        compare_reports(second, first)
    second["pmf"] = None
    assert compare_reports(second, first)["passed"] is None


def test_figures_preserve_gaps_and_readiness(tmp_path: Path) -> None:
    data = report()
    plot_report(data, tmp_path / "pmf.png")
    assert (tmp_path / "pmf.png").stat().st_size > 1000
    with pytest.raises(FileExistsError):
        plot_report(data, tmp_path / "pmf.png")
    data.pop("overlap")
    data["pmf"] = None
    plot_report(data, tmp_path / "readiness.png")
    assert (tmp_path / "readiness.png").stat().st_size > 1000


def test_uncensored_chemical_review_retains_owner_swaps() -> None:
    blocks = [
        dict(
            window="w",
            values=np.array([[0, 0, 0, 1.0], [0, 0, 0, 2.0]]),
            original_rows=[dict(owner_changed="true"), dict(owner_changed="false")],
        )
    ]
    result = chemical_diagnostics(blocks, ["cn"])
    assert result["samples_censored"] == 0 and result["label_swap_present"]
    assert result["split_block_diagnostics"][0]["mean_shift"] == 1.0


def test_contiguous_bootstrap_is_reproducible_and_keeps_failed_support(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rng = np.random.default_rng(19)
    values = np.column_stack([np.cumsum(rng.normal(size=200)), np.zeros((200, 2))])
    block = dict(window="w", values=values, native_bias_kcal_mol=np.zeros(200))
    windows = [
        dict(
            label="w",
            bias=dict(
                feature_names=("q_acyl", "q_act"),
                centers_A=(0.0, 0.0),
                force_constants=(1.0, 1.0),
                outer_extent_A=10.0,
            ),
        )
    ]

    def estimator(windows: list, blocks: list, **kwargs: object) -> dict:
        mean = float(np.mean(blocks[0]["values"]))
        return dict(
            blocks=[dict(discarded_frames=5, statistical_inefficiency=3.0)],
            pmf=dict(
                minimum_to_barrier_kcal_mol=mean,
                basin_minimum_difference_kcal_mol=mean / 2,
            ),
        )

    monkeypatch.setattr(bootstrap_module, "analyze", estimator)
    options = dict(analysis_options={}, replicates=20, chunk_frames=10, seed=12)
    first = bootstrap_module.bootstrap(windows, [block], **options)
    second = bootstrap_module.bootstrap(windows, [block], **options)
    assert first == second and len(first["estimates"]) == 20
    assert first["halfwidth_kcal_mol"][0] > 0
    with pytest.raises(ValueError, match="twenty"):
        bootstrap_module.bootstrap(windows, [block], **(options | dict(replicates=2)))
    with pytest.raises(ValueError, match="shorter"):
        bootstrap_module.bootstrap(windows, [block], **(options | dict(chunk_frames=2)))
    with pytest.raises(ValueError, match="exceed"):
        bootstrap_module.bootstrap(
            windows, [block], **(options | dict(chunk_frames=1000))
        )
    monkeypatch.setattr(
        bootstrap_module,
        "analyze",
        lambda *a, **kw: dict(
            blocks=[dict(discarded_frames=0, statistical_inefficiency=1.0)], pmf=None
        ),
    )
    failed = bootstrap_module.bootstrap(windows, [block], **options)
    assert failed["failed_support_replicates"] == 20 and not failed["passed"]


def test_local_watcher_detects_all_window_coverage(tmp_path: Path) -> None:
    import json

    script = (
        Path(__file__).resolve().parents[1] / "scripts/watch_relay_umbrella_analysis.py"
    )
    spec = importlib.util.spec_from_file_location("relay_watcher", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cfg = dict(windows=[dict(label="a"), dict(label="b")])
    assert module.discover(cfg, tmp_path, "analysis_acceptance.json") == ([], False)
    for label in ("a", "b"):
        folder = tmp_path / label
        folder.mkdir()
        (folder / "analysis_acceptance.json").write_text(
            json.dumps(dict(window=label, accepted=True, sampled_duration_ps=5))
        )
    assert module.discover(cfg, tmp_path, "analysis_acceptance.json")[1]
    with pytest.raises(ValueError, match="Unexpected"):
        module.discover(
            dict(windows=[dict(label="a")]), tmp_path, "analysis_acceptance.json"
        )


def test_watcher_once_writes_readiness_and_completion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import argparse
    import json

    script = (
        Path(__file__).resolve().parents[1] / "scripts/watch_relay_umbrella_analysis.py"
    )
    spec = importlib.util.spec_from_file_location("relay_watcher_once", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    request = tmp_path / "request.json"
    request.write_text(
        json.dumps(dict(windows=[dict(label="a")], accepted_block_receipts=[]))
    )
    receipts = tmp_path / "receipts"
    receipts.mkdir()
    events = []
    monkeypatch.setattr(module, "journal", lambda *args, **kwargs: events.append(args))
    monkeypatch.setattr(
        module,
        "analyze_request",
        lambda path: dict(status="readiness_no_accepted_data"),
    )
    monkeypatch.setattr(
        module, "plot_report", lambda report, path: path.write_text("plot stub")
    )
    args = argparse.Namespace(
        request=request,
        output_root=tmp_path / "watcher",
        receipts_root=receipts,
        receipt_pattern="analysis_acceptance.json",
        max_minutes=1,
        poll_seconds=60,
        once=True,
        resume_output=False,
        journal_root=tmp_path,
    )
    module.watch(args)
    completion = json.loads((args.output_root / "watcher_completion.json").read_text())
    assert completion["reason"] == "once" and not completion["full_coverage_analyzed"]
    assert (tmp_path / "readiness_update_001.json").is_file()
    assert events[0][0] == "intent" and events[-1][0] == "outcome"
