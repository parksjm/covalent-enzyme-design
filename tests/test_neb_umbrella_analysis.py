"""Completed-block gates, uncensored statistics, unbiased MBAR and missing-bin plots."""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.neb_umbrella_analysis import (
    analyze,
    collect,
    components,
    digest,
    histogram_fes,
)
from scripts.analyze_neb_umbrella import refresh


def campaign(tmp_path: Path, count: int = 2) -> Path:
    """Create known synthetic biased windows and production completion receipts."""
    rng = np.random.default_rng(4)
    windows = []
    for i in range(count):
        label = f"window_{i:03d}"
        stage = "production_01"
        bias = dict(
            q_acyl_center=(-0.3 if i == 0 else 0.3),
            q_h_center=0.0,
            k_acyl=2.0,
            k_h=2.0,
            outer_extent=10.0,
        )
        windows.append(
            dict(
                label=label,
                bias=bias,
                stages=[
                    dict(label="heat", kind="heat"),
                    dict(label=stage, kind="production"),
                    dict(label="partial", kind="production"),
                ],
            )
        )
        directory = tmp_path / label / stage
        folder = directory / "observables"
        folder.mkdir(parents=True)
        data = rng.normal(size=(80, 4)) * 0.1
        data[:, 0] += bias["q_acyl_center"]
        rows = []
        for frame, s in enumerate(data):
            rows.append(
                dict(
                    frame=frame,
                    q_acyl=s[0],
                    q_h=s[1],
                    cn_progress=s[2],
                    cn_nalpha_hpool=s[3] + 2,
                    hg1_owner=13916 if frame < 40 else 19929,
                    nalpha_h1_owner=13905,
                    nalpha_h2_owner=13905,
                    owner_changed=bool(frame % 2),
                )
            )
        p = folder / "observables.csv"
        with p.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (folder / "observables.json").write_text(
            json.dumps(dict(csv_sha256=digest(p), frames=80, frames_censored=0))
        )
        (directory / "acceptance.json").write_text(
            json.dumps(dict(finite_energy_temperature=True, completed_steps=10000))
        )
        np.savez(directory / "samples_and_cross_bias.npz", cv_values=data[:, :2])
    prepared = tmp_path / "prepared.json"
    prepared.write_text(json.dumps(dict(temperature_K=300, windows=windows)))
    (tmp_path / "prepared.sha256").write_text(digest(prepared))
    return prepared


def test_collection_gates_only_complete_production_and_hashes(tmp_path: Path) -> None:
    p = campaign(tmp_path)
    cfg, windows = collect(p)
    assert len(windows) == 2 and len(windows[0]["blocks"]) == 1
    csv_path = tmp_path / "window_000/production_01/observables/observables.csv"
    csv_path.write_text(csv_path.read_text() + "\n")
    with pytest.raises(ValueError, match="identity"):
        collect(p)
    p.write_text("{}")
    with pytest.raises(ValueError, match="identity"):
        collect(p)


def test_exact_mbar_diagnostics_preserve_branches_and_missing_bins(
    tmp_path: Path,
) -> None:
    cfg, windows = collect(campaign(tmp_path))
    result = analyze(cfg, windows, bins=8)
    assert result["windows_sampled"] == 2 and result["chemical_frames_censored"] == 0
    assert result["window_summaries"][0]["changed_owner_frames"] == 40
    assert len(result["window_summaries"][0]["owner_states"]) == 2
    assert result["mbar_status"] == "interim_bias_only_common_U0"
    assert 0 < result["unbiased_effective_sample_size"] <= sum(result["mbar_N_k"])
    assert result["primary_fes"]["unsampled_bins"] > 0
    assert sum(sum(r) for r in result["primary_fes"]["probability"]) == pytest.approx(
        1.0
    )
    assert result["activation_barrier_claim"] is False
    cfg["windows"].append(dict(label="missing", bias=cfg["windows"][0]["bias"]))
    partial = analyze(cfg, windows, bins=4)
    assert partial["mbar_N_k"][-1] == 0 and partial["status"] == "interim_incomplete"
    assert analyze(cfg, [], bins=4)["status"] == "no_completed_production"
    assert analyze(cfg, windows[:1], bins=4)["mbar_status"].startswith("insufficient")
    with pytest.raises(ValueError):
        analyze(cfg, windows, bins=1)


def test_refresh_plots_and_unchanged_input_skip(tmp_path: Path) -> None:
    p = campaign(tmp_path)
    out = tmp_path / "analysis"
    first = refresh(p, out, bins=5)
    second = refresh(p, out, bins=5)
    assert first["generation"] == second["generation"]
    assert (out / "latest/analysis.json").is_file()
    assert len(list(out.glob("update_*"))) == 1
    assert components(np.eye(3)) == [[0], [1], [2]]
    fes = histogram_fes(
        np.array([[0.0, 0.0], [1.0, 1.0]]), np.array([0.5, 0.5]), 300, 2
    )
    assert fes["unsampled_bins"] == 2


def test_cli_watch_hook_and_invalid_walltime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The bounded watcher invokes analysis only and rejects invalid walltime."""
    from scripts import analyze_neb_umbrella as cli
    import sys

    calls = []
    monkeypatch.setattr(cli, "refresh", lambda *args: calls.append(args))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "analyze",
            "--prepared",
            str(tmp_path / "prepared.json"),
            "--watch",
            "--hours",
            "0.001",
        ],
    )
    clock = iter([0.0, 10.0])
    monkeypatch.setattr(cli.time, "monotonic", lambda: next(clock))
    cli.main()
    assert len(calls) == 1
    monkeypatch.setattr(
        sys, "argv", ["analyze", "--prepared", "prepared.json", "--hours", "-1"]
    )
    with pytest.raises(SystemExit):
        cli.main()


def test_disconnected_sampled_overlap_suppresses_global_fes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A disconnected MBAR graph must not publish either globally normalized FES."""
    import pymbar

    cfg, windows = collect(campaign(tmp_path))
    original = pymbar.MBAR.compute_overlap
    monkeypatch.setattr(
        pymbar.MBAR, "compute_overlap", lambda self: {"matrix": np.eye(len(self.N_k))}
    )
    result = analyze(cfg, windows, bins=5)
    assert result["status"] == "interim_disconnected_overlap"
    assert result["fes_status"] == "suppressed_disconnected_sampled_overlap"
    assert result["connected"] is False
    assert "primary_fes" not in result and "diagnostic_fes" not in result
    assert "coverage_points" in result and "mbar_overlap" in result
    assert result["relative_basin_claim"] is False
    monkeypatch.setattr(pymbar.MBAR, "compute_overlap", original)
