"""Synthetic validation of bounded targeted analysis decisions."""

import csv
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
stats = importlib.import_module("relay_targeted_statistics")
watch = importlib.import_module("watch_targeted_relay_analysis")


def write_csv(path: Path, times: list[float]) -> None:
    path.with_name("run.mdin").write_text("dt=.0005,ntwx=1000,nstlim=10000")
    with path.open("w") as h:
        writer = csv.DictWriter(
            h,
            fieldnames=[
                "frame",
                "time_ps",
                *stats.FEATURES,
                *stats.OWNERS,
                "hydrogen_owners",
            ],
        )
        writer.writeheader()
        for i, t in enumerate(times):
            row = dict(frame=i + 1, time_ps=t, hydrogen_owners="all retained")
            row.update({k: (i % 3) * 0.1 for k in stats.FEATURES})
            row.update({k: "Nalpha" for k in stats.OWNERS})
            writer.writerow(row)


def test_discard_once(tmp_path: Path) -> None:
    paths = [tmp_path / f"{i}.csv" for i in range(2)]
    for i, p in enumerate(paths):
        write_csv(p, [i * 5 + j * 0.5 for j in range(1, 11)])
    full, kept = stats.production_rows(paths)
    assert len(full) == 20 and len(kept) == 10
    result = stats.summarize(kept)
    assert result["frames"] == 10
    assert result["owners"]["hg1_owner"] == {"Nalpha": 10}
    assert all(v["ess"] > 0 for v in result["features"].values())
    constant = [dict(r, **{k: 1 for k in stats.FEATURES}) for r in kept]
    assert stats.summarize(constant)["features"]["q_del_A"]["ess"] is None
    with pytest.raises(ValueError):
        stats.summarize(kept[:2])
    kept[0]["q_del_A"] = "nan"
    with pytest.raises(ValueError):
        stats.summarize(kept)
    with pytest.raises(ValueError):
        stats.production_rows(paths[::-1])
    with pytest.raises(ValueError):
        stats.production_rows([])


def test_overlap_and_bias() -> None:
    a = np.array([-0.2, 0.1, 0.2, 0.3])
    assert stats.histogram_overlap(a, a) == pytest.approx(1)
    assert stats.histogram_overlap(a, a + 5) == 0
    assert stats.histogram_overlap(np.c_[a, a], np.c_[a, a]) == pytest.approx(1)
    stats.require_same_bias({"centers": [1, 2]}, {"centers": [1, 2]})
    with pytest.raises(ValueError):
        stats.require_same_bias({"centers": [1, 2]}, {"centers": [1, 3]})


@pytest.mark.parametrize(
    "kwargs,prefix",
    [
        ({"chemistry_failure": True}, "inspect_chemistry"),
        ({"ess": None}, "longer_sampling"),
        ({"ess": 10}, "longer_sampling"),
        ({"drift": 0.3}, "longer_sampling"),
        ({"same_bias": True, "biased_overlap": 0.8, "overlap": 0.3}, "orthogonal"),
        ({"overlap": 0.1}, "additional_centers"),
        ({}, "conditional"),
    ],
)
def test_recommendations(kwargs: dict, prefix: str) -> None:
    args = dict(
        ess=100, drift=0.01, overlap=0.8, same_bias=False, biased_overlap=0.8, gap=0.01
    )
    args.update(kwargs)
    assert stats.recommendation(**args).startswith(prefix)


def test_native_success_outer_failure(tmp_path: Path) -> None:
    assert not watch.stage_ready(tmp_path, 100)
    receipts = {
        "acceptance.json": dict(
            accepted=True,
            completed_steps=100,
            finite_energy_temperature=True,
            every_integration_step_SCC_converged=True,
            checkpoint_sha256="abc",
        ),
        "process.json": dict(exit_status=0, ended_utc="now"),
        "analysis_outcome.json": dict(status="failed", error="derived only"),
        "wrapper.json": dict(exit_status=1, error="unsupported pidfd"),
    }
    for n, v in receipts.items():
        (tmp_path / n).write_text(json.dumps(v))
    assert watch.stage_ready(tmp_path, 100)
    assert not watch.stage_ready(tmp_path, 200)
    (tmp_path / "process.json").write_text(
        json.dumps(dict(exit_status=1, ended_utc="now"))
    )
    assert not watch.stage_ready(tmp_path, 100)


def test_complete_report(tmp_path: Path) -> None:
    reporting = importlib.import_module("relay_targeted_report")
    labels = [
        "activation_021_from_022",
        "activation_022_from_021",
        "delivery_027_qdel_1p50",
        "delivery_028_qdel_1p75",
    ]
    parents = ["activation_full_021", "activation_full_022", "delivery_full_026"]
    bias = dict(
        feature_names=["q_acyl", "q_act"], centers_A=[1, 2], force_constants=[25, 25]
    )
    (tmp_path / "prepared.json").write_text(
        json.dumps(dict(windows=[dict(label=p, bias=bias) for p in parents]))
    )
    rows = []
    for i, label in enumerate(labels):
        stages = []
        for j in range(3, 5):
            p = (
                tmp_path
                / "additional_sampling_20261126"
                / label
                / f"{j:02}_production_5ps"
                / "raw_features.csv"
            )
            p.parent.mkdir(parents=True)
            write_csv(p, [5 * (j - 3) + 0.5 * k for k in range(1, 11)])
            mdin = p.with_name("run.mdin")
            mdin.write_text("irest=0,ntx=1,dt=.0005,ntwx=1000,nstlim=10000")
            stages.append(
                dict(
                    kind="production",
                    input=str(
                        mdin.relative_to(tmp_path / "additional_sampling_20261126")
                    ),
                    seed=j,
                )
            )
        rows.append(
            dict(
                label=label,
                bias=bias,
                stages=stages,
                source_lineage=dict(opposite_start_target=parents[i % 2]),
            )
        )
    for label in parents:
        for j in range(3, 9):
            p = tmp_path / label / f"{j:02}_production_5ps" / "raw_features.csv"
            p.parent.mkdir(parents=True)
            write_csv(p, [5 * (j - 3) + 0.5 * k for k in range(1, 11)])
    report = reporting.make_report(tmp_path, dict(windows=rows), tmp_path)
    assert len(report["comparisons"]) == 4
    assert report["windows"][labels[0]]["summary"]["frames"] == 10
    assert report["windows"][parents[0]]["summary"]["frames"] == 50
    assert (tmp_path / "targeted_relay_activation_relaxation.png").stat().st_size > 1000


def test_watcher_success_and_dead_workers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import hashlib
    import parmed
    import nyl12_covalent.relay_umbrella_runtime as runtime

    reporting = importlib.import_module("relay_targeted_report")
    prepared = tmp_path / "prepared.json"
    cfg = dict(
        frozen_files_sha256={},
        topology="unused",
        windows=[dict(label="one", stages=[dict(input="one/stage/run.mdin", steps=1)])],
    )
    prepared.write_text(json.dumps(cfg))
    prepared.with_suffix(".sha256").write_text(
        hashlib.sha256(prepared.read_bytes()).hexdigest()
    )
    monkeypatch.setattr(watch, "journal", lambda *args: None)
    monkeypatch.setattr(watch, "stage_ready", lambda *args: True)
    monkeypatch.setattr(parmed, "load_file", lambda *args: object())
    monkeypatch.setattr(runtime, "ensure_diagnostics", lambda *args: True)
    monkeypatch.setattr(reporting, "make_report", lambda *args: {})
    watch.run(prepared, tmp_path, tmp_path, 1, 10)
    monkeypatch.setattr(watch, "stage_ready", lambda *args: False)
    monkeypatch.setattr(watch, "workers_alive", lambda *args: False)
    with pytest.raises(RuntimeError, match="no live"):
        watch.run(prepared, tmp_path, tmp_path, 1, 10)
    assert (
        json.loads((tmp_path / "targeted_relay_partial.json").read_text())["status"]
        == "partial_or_failed"
    )
    with pytest.raises(ValueError):
        watch.run(prepared, tmp_path, tmp_path, 25, 10)


def test_process_scan_journal_and_cli(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert not watch.workers_alive(tmp_path / "nonexistent_campaign")
    calls = []
    monkeypatch.setattr(watch.subprocess, "run", lambda *a, **k: calls.append(a))
    watch.journal(tmp_path, "intent", "test synthetic journal command")
    assert calls
    monkeypatch.setattr(watch, "run", lambda *a: calls.append(a))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "watch",
            "--prepared",
            str(tmp_path / "prepared.json"),
            "--output",
            str(tmp_path),
            "--journal-root",
            str(tmp_path),
        ],
    )
    watch.main()
    assert len(calls) == 2


def test_interstage_wrapper_identity_without_sander(tmp_path: Path) -> None:
    import os

    proc = Path("/proc") / str(os.getpid())
    ticks = int((proc / "stat").read_text().rsplit(")", 1)[1].split()[19])
    status = tmp_path / "one/wrapper_state/wrapper_status.json"
    status.parent.mkdir(parents=True)
    status.write_text(
        json.dumps(dict(pid=os.getpid(), start_ticks=ticks, status="starting"))
    )
    assert watch.workers_alive(tmp_path)
    status.write_text(
        json.dumps(dict(pid=os.getpid(), start_ticks=ticks - 1, status="starting"))
    )
    assert not watch.workers_alive(tmp_path)


def test_native_float32_exact_burn_boundary(tmp_path: Path) -> None:
    paths = [tmp_path / str(i) / "raw_features.csv" for i in range(2)]
    for i, p in enumerate(paths):
        p.parent.mkdir()
        write_csv(p, [float(np.float32(6 + i * 5 + j * 0.05)) for j in range(1, 101)])
        p.with_name("run.mdin").write_text("dt=.0005,ntwx=100,nstlim=10000")
    full, kept = stats.production_rows(paths)
    assert len(full) == 200 and len(kept) == 100
    assert full[99]["production_time_ps"] == 5.0
    assert kept[0]["production_time_ps"] == 5.05
    rows = paths[1].read_text().replace(str(float(np.float32(11.05))), "11.06", 1)
    paths[1].write_text(rows)
    with pytest.raises(ValueError, match="clock discontinuity"):
        stats.production_rows(paths)
