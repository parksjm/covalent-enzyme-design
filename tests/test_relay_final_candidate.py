"""Actual final-horizon gates and count-preserving uncertainty regressions."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent import relay_final_bootstrap as boot
from nyl12_covalent import relay_final_sources as source
from nyl12_covalent import relay_final_statistics as stats
from nyl12_covalent import relay_final_analysis as analysis


def test_circular_arbitrary_lengths_and_tail() -> None:
    counts = np.array([1100, 500])
    seen = set()
    for seed in range(15):
        ids = boot.circular_indices(counts, 200, seed)
        assert len(ids) == 1600
        assert ids[:1100].min() >= 0 and ids[:1100].max() < 1100
        assert ids[1100:].min() >= 1100 and ids[1100:].max() < 1600
        seen.update(ids.tolist())
    assert set(range(1600)) == seen
    assert np.array_equal(
        boot.circular_indices(counts, 100, 1), boot.circular_indices(counts, 100, 1)
    )
    for bad, width in [(np.array([100]), 100), (counts, 50), (np.array([[500]]), 100)]:
        with pytest.raises(ValueError):
            boot.circular_indices(bad, width, 1)


def test_paired_circular_fit_and_disconnection(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(boot, "_WORK", None)
    with pytest.raises(RuntimeError):
        boot.draw((100, 1))
    calls = []

    def fit(u: np.ndarray, n: np.ndarray) -> tuple:
        calls.append(u[0].copy())
        return np.zeros(u.shape[1]), {"connected": len(calls) % 2 == 0}

    def evaluate(lw: np.ndarray, p: dict, ids: np.ndarray, support_stats: bool) -> dict:
        return {name: {"pmf": np.array([float(ids[0]), 0.0])} for name in p}

    monkeypatch.setattr(boot, "fit_logweights", fit)
    monkeypatch.setattr(boot, "evaluate", evaluate)
    result = boot.circular_bootstrap(
        np.arange(500)[None, :], np.array([500]), {"a": {}, "b": {}}, 200, workers=1
    )
    assert len(calls) == 200 and result["last5ps_dropped"] is False
    assert result["projections"]["a"] == result["projections"]["b"]
    assert result["projections"]["a"]["lower_95"] == [None, None]
    monkeypatch.setattr(
        boot, "fit_logweights", lambda u, n: (np.zeros(500), {"connected": True})
    )
    result = boot.circular_bootstrap(
        np.arange(500)[None, :], np.array([500]), {"a": {}}, 100, workers=1
    )
    assert result["projections"]["a"]["lower_95"][1] == 0
    with pytest.raises(ValueError):
        boot.circular_bootstrap(np.zeros((1, 500)), np.array([500]), {}, 100, workers=5)
    monkeypatch.setattr(boot, "paired_bootstrap", lambda *args, **kw: {"same": True})
    assert boot.nonoverlap_baseline(np.zeros((1, 500)), np.array([500]), {})[
        "projections"
    ]["same"]


def test_actual_source_gate_and_no45_substitution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared = tmp_path / "prepared.json"
    prepared.write_text(json.dumps({"windows": [{"label": "w"}]}))
    monkeypatch.setattr(source, "COHORTS", ((prepared, 10, "test", 1),))
    assert not source.readiness()["ready"] and source.final_data() is None
    for stage in (3, 4):
        p = tmp_path / "w" / f"{stage:02d}_production_5ps/analysis_acceptance.json"
        p.parent.mkdir(parents=True)
        p.write_text(json.dumps({"accepted": True, "window": "w"}))
    assert source.readiness()["ready"]
    monkeypatch.setattr(
        source,
        "load_cohort",
        lambda *args: dict(states=[], counts=np.array([]), proofs=[], samples=[]),
    )
    with pytest.raises(ValueError, match="actualmain60"):
        source.final_data()
    p.write_text(json.dumps({"accepted": False, "window": "w"}))
    with pytest.raises(ValueError):
        source.readiness()


def test_frozen_physical_grids() -> None:
    cfg = json.loads(
        (
            source.FIXED
            / "automatic_progression_20261126/checkpoint_45ps/analysis_request.json"
        ).read_text()
    )
    q = np.asarray(cfg["path_3cv_A"])
    grids = analysis.definitions(dict(q=q, config=cfg))
    assert grids["delivery_2D"]["shape"] == (90, 60)
    assert "BGP_lambda_x0.5" in grids and "BGP_lambda_x2.0" in grids
    assert stats.intersection(np.ones(100), np.ones(100)) == 1
    assert stats.intersection(np.zeros(100), np.ones(100)) == 0


def test_six_derived_sync_has_no_native_files() -> None:
    spec = importlib.util.spec_from_file_location(
        "finalwatch", Path("scripts/watch_relay_final_candidate.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    binding = dict(
        host="azalea",
        reviewed=True,
        prepared_sha256=module.digest(module.MIXING),
        ssh_key="/home/jpc/.ssh/id_ed25519_aster",
        remote_root=str(module.MIXING.parent),
    )
    commands = module.sync_commands(binding)
    assert len(commands) == 6
    assert all("--exclude=*" in c and "--include=*.csv" in c for c in commands)
    assert all("azalea:" in c[-2] and ".nc" not in " ".join(c) for c in commands)


def test_physical_independent_start_profiles(tmp_path: Path) -> None:
    import csv

    row = dict(
        frame=1,
        r_c_n=1.5,
        r_c_og1=2.5,
        r_og1_hg1=1.0,
        q_act_A=-1.0,
        r_nleave_hg1=3.0,
        hg1_owner="13916",
    )
    csvfile = tmp_path / "features.csv"
    with csvfile.open("w") as stream:
        w = csv.DictWriter(stream, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    receipt = tmp_path / "analysis_acceptance.json"
    receipt.write_text(json.dumps({"csv_file": csvfile.name}))
    labels = [
        "main60/activation_full_026",
        "main60/activation_full_028",
        "main60/delivery_full_008",
        "main60/delivery_full_009",
    ]
    labels += [
        "mixing_checks30/" + x
        for x in [
            "A026_original",
            "A026_late",
            "A028_original",
            "A028_late",
            "D008_from_D010",
            "D009_from_D007",
        ]
    ]
    data = dict(
        states=[dict(state_id=x, window_statistics={}) for x in labels],
        counts=np.ones(10, dtype=int),
        samples=[dict(receipt=str(receipt), source_frame1based=1)] * 10,
    )
    result = stats.summarize(data)
    assert len(result["independent_start_pairs"]) == 6
    s = result["lineages"][labels[0]]
    assert s["distances"]["Nalpha_HG1"]["mean"] == 2.0
    assert s["hg1_owner_fractions"] == {"13916": 1.0}
    data["states"][0]["state_id"] = "wrong"
    with pytest.raises(ValueError, match="lineage missing"):
        stats.summarize(data)


def test_full_final_point_outputs_without_extra_fit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = json.loads(
        (
            source.FIXED
            / "automatic_progression_20261126/checkpoint_45ps/analysis_request.json"
        ).read_text()
    )
    q = np.asarray(cfg["path_3cv_A"])
    data = dict(q=q, counts=np.array([len(q)]), config=cfg, states=[], samples=[])

    def fake_point(d: dict) -> tuple:
        return (
            dict(states=100, frames=len(q), counts=[len(q)]),
            {},
            np.zeros((1, len(q))),
            np.full(len(q), -np.log(len(q))),
        )

    monkeypatch.setattr(analysis, "point", fake_point)
    monkeypatch.setattr(
        analysis,
        "fit_logweights",
        lambda u, n: (np.full(u.shape[1], -np.log(u.shape[1])), {"connected": True}),
    )
    monkeypatch.setattr(analysis, "summarize", lambda d: {"checked": True})
    calls = []
    monkeypatch.setattr(
        analysis,
        "circular_bootstrap",
        lambda u, n, p, b, w: calls.append(b) or {"block": b},
    )
    monkeypatch.setattr(
        analysis, "nonoverlap_baseline", lambda *args: {"baseline": True}
    )
    out = tmp_path / "candidate"
    acceptance = analysis.analyze(data, out)
    assert calls == [100, 200] and acceptance["primary_replacement"] is False
    assert (out / "final_candidate_point.png").is_file()
    assert (out / "circular_10ps_conditional95.json").is_file()
    assert (out / "matched_final_weights.npz").is_file()


def test_final_source_full100_exact_counts(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(source, "readiness", lambda: {"ready": True})

    def load(path: Path, horizon: int, reference: dict, group: str, n: int) -> dict:
        return dict(
            states=[dict(state_id=f"{group}/{i}") for i in range(n)],
            counts=np.full(n, (horizon - 5) * 20),
            q=np.zeros((n * (horizon - 5) * 20, 3)),
            proofs=[dict(source_sha256={f"{group}.csv": group})],
            samples=[],
            prepared_path=str(path),
            prepared_sha256="proof",
        )

    monkeypatch.setattr(source, "load_cohort", load)
    data = source.final_data()
    assert data["q"].shape == (84200, 3) and len(data["states"]) == 100
    assert data["main_actual_horizon_ps"] == 60


def test_frozen_watcher_gates_then_runs_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = importlib.util.spec_from_file_location(
        "watchfinal", Path("scripts/watch_relay_final_candidate.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    prepared = tmp_path / "prepared.json"
    prepared.write_text("{}")
    module.COHORTS = ((prepared, 60, "main", 57),)
    manifest = tmp_path / "source_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "files": {},
                "static_inputs": {},
                "external_inputs": {str(prepared): module.digest(prepared)},
            }
        )
    )
    monkeypatch.setattr(
        module, "sync_commands", lambda binding: [["derived-only-sync"]]
    )
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: None)
    states = iter([{"ready": False}, {"ready": True}])
    monkeypatch.setattr(module, "readiness", lambda: next(states))
    monkeypatch.setattr(module, "final_data", lambda: {"loaded_actual60": True})
    monkeypatch.setattr(module.time, "sleep", lambda *args: None)
    binding = tmp_path / "host_binding.json"
    binding.write_text("{}")
    calls = []
    monkeypatch.setattr(
        module,
        "analyze",
        lambda data, path, workers: calls.append((data, path, workers)),
    )
    module.watch(tmp_path, manifest, binding)
    assert len(calls) == 1 and calls[0][2] == 4
    status = json.loads((tmp_path / "monitor_state.json").read_text())
    assert (
        status["status"] == "candidate_complete"
        and status["primary_replacement"] is False
    )
    with pytest.raises(ValueError):
        module.watch(tmp_path, manifest, binding, hours=100)


def test_watcher_cli_and_binding_rejections(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import sys

    spec = importlib.util.spec_from_file_location(
        "finalcli", Path("scripts/watch_relay_final_candidate.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with pytest.raises(ValueError, match="reviewed"):
        module.sync_commands(dict(host="unapproved", reviewed=False))
    with pytest.raises(ValueError, match="identity"):
        module.sync_commands(
            dict(host="geranium", reviewed=True, prepared_sha256="bad")
        )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "watch",
            "--output",
            str(tmp_path),
            "--source-manifest",
            str(tmp_path / "manifest.json"),
            "--host-binding",
            str(tmp_path / "binding.json"),
        ],
    )
    calls = []
    monkeypatch.setattr(module, "watch", lambda *args: calls.append(args))
    module.main()
    assert len(calls) == 1 and calls[0][-2:] == (72, 120)
