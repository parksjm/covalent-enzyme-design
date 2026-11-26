"""Acceptance and provenance failures stop new relay execution explicitly."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from nyl12_covalent.relay_umbrella_runtime import native_acceptance, verify_plan
from nyl12_covalent.neb_umbrella_campaign import digest


def native_text(steps: int = 2) -> str:
    """Return minimal step-resolved native output."""
    return "\n".join(
        [
            *(
                f"QMMM SCC-DFTB: SCC-DFTB for step {n} converged in 20 cycles."
                for n in range(steps + 1)
            ),
            f"NSTEP = {steps} TIME(PS) = 0.001 TEMP(K) = 300.0",
            "EPtot = -12.0",
        ]
    )


def test_acceptance_requires_resolved_scc(tmp_path: Path) -> None:
    (tmp_path / "run.mdout").write_text(native_text())
    assert native_acceptance(tmp_path, 2, [250, 350])["completed_steps"] == 2
    (tmp_path / "run.mdout").write_text(
        native_text().replace("step 1 converged", "step 1 unresolved")
    )
    with pytest.raises(ValueError, match="Unresolved SCC"):
        native_acceptance(tmp_path, 2, [250, 350])


@pytest.mark.parametrize(
    "replacement,match",
    [
        ("NSTEP = 1", "expected step"),
        ("TEMP(K) = 600.0", "temperature outside"),
        ("EPtot = NaN", "electronic history"),
    ],
)
def test_acceptance_rejects_bad_native(
    tmp_path: Path, replacement: str, match: str
) -> None:
    text = native_text()
    target = (
        "NSTEP = 2"
        if replacement.startswith("NSTEP")
        else "TEMP(K) = 300.0"
        if replacement.startswith("TEMP")
        else "EPtot = -12.0"
    )
    (tmp_path / "run.mdout").write_text(text.replace(target, replacement))
    with pytest.raises(ValueError, match=match):
        native_acceptance(tmp_path, 2, [250, 350])


def test_plan_pins_inputs_runtime_and_source(tmp_path: Path) -> None:
    data = tmp_path / "input"
    data.write_text("unchanged")
    prepared = tmp_path / "prepared.json"
    prepared.write_text(json.dumps({"frozen_files_sha256": {"input": digest(data)}}))
    plan = tmp_path / "plan.json"
    plan.write_text(
        json.dumps(
            dict(
                prepared=str(prepared),
                prepared_sha256=digest(prepared),
                runtime_sha256={str(data): digest(data)},
                source_sha256={str(data): digest(data)},
            )
        )
    )
    plan.with_suffix(".sha256").write_text(digest(plan))
    assert verify_plan(plan)[2] == tmp_path
    data.write_text("changed")
    with pytest.raises(ValueError, match="Frozen scientific input changed"):
        verify_plan(plan)


@pytest.mark.parametrize("corruption", ["checkpoint", "lineage"])
def test_segment_runs_and_revalidates_checkpoint_prefix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, corruption: str
) -> None:
    """A second block revalidates, never reruns, its accepted predecessor."""
    import sys
    from types import SimpleNamespace
    import numpy as np
    from nyl12_covalent import relay_umbrella_runtime as runtime

    stages = []
    for index, kind in enumerate(("heat", "equil", "production", "production")):
        label = f"stage{index}"
        (tmp_path / "window" / label).mkdir(parents=True)
        (tmp_path / "window" / label / "run.mdin").write_text("input")
        (tmp_path / "window" / label / "umbrella.RST").write_text("bias")
        stages.append(
            dict(
                label=label,
                kind=kind,
                input=f"window/{label}/run.mdin",
                steps=2,
                seed=index,
            )
        )
    (tmp_path / "inputs").mkdir()
    (tmp_path / "inputs/system.prmtop").write_text("top")
    (tmp_path / "seed").write_text("seed")
    plan_path = tmp_path / "plan.json"
    plan_path.write_text("plan")
    cfg = dict(
        windows=[dict(label="window", seed="seed", stages=stages)],
        system=dict(active_mask="mobile"),
    )
    plan = dict(
        sander="/sander",
        environment={},
        heating_temperature_gate_K=[0, 400],
        temperature_gate_K=[250, 350],
    )
    monkeypatch.setattr(runtime, "verify_plan", lambda path: (plan, cfg, tmp_path))
    monkeypatch.setattr(runtime.parmed, "load_file", lambda path: object())
    initial = SimpleNamespace(
        coordinates=np.zeros((104922, 3)), box=np.array([1, 1, 1, 90, 90, 90])
    )
    monkeypatch.setattr(runtime.Rst7, "open", lambda path: initial)
    monkeypatch.setattr(
        runtime,
        "AmberMask",
        lambda top, mask: SimpleNamespace(Selected=lambda: range(4165)),
    )
    monkeypatch.setattr(
        runtime,
        "_checkpoint",
        lambda path, *args: ({"checkpoint_sha256": digest(path)}, np.zeros((1, 3))),
    )
    calls = []

    def run(stage: Path, *args: object) -> Path:
        calls.append(stage.name)
        (stage / "process.json").write_text(json.dumps(dict(exit_status=0)))
        (stage / "run.mdout").write_text(native_text())
        (stage / "trajectory.nc").write_text("native coordinates")
        (stage / "restart.ncrst").write_text(stage.name)
        return stage / "restart.ncrst"

    monkeypatch.setattr(runtime, "_run_stage", run)

    def diagnostic(cfg: dict, row: dict, stage: Path, top: object) -> None:
        (stage / "analysis_request.json").write_text("{}")
        (stage / "raw_features.csv").write_text("features")
        (stage / "samples_and_cross_bias.npz").write_text("matrix")

    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.relay_umbrella_diagnostics",
        SimpleNamespace(block_diagnostics=diagnostic),
    )
    runtime.run_segment(plan_path, "window", 1)
    assert calls == ["stage0", "stage1", "stage2"]
    runtime.run_segment(plan_path, "window", 2)
    assert calls == ["stage0", "stage1", "stage2", "stage3"]
    assert (
        json.loads((tmp_path / "window/segment_receipt_limit_02.json").read_text())[
            "production_blocks"
        ]
        == 2
    )
    if corruption == "checkpoint":
        (tmp_path / "window/stage0/restart.ncrst").write_text("corrupted")
        match = "checkpoint changed"
    else:
        path = tmp_path / "window/stage3/input_identity.json"
        identity = json.loads(path.read_text())
        identity["sha256"][str(tmp_path / "window/stage2/restart.ncrst")] = (
            "wrong-parent"
        )
        path.write_text(json.dumps(identity))
        match = "predecessor lineage changed"
    with pytest.raises(ValueError, match=match):
        runtime.run_segment(plan_path, "window", 2)
    assert (
        json.loads((tmp_path / "window/window_status.json").read_text())["status"]
        == "failed"
    )
    with pytest.raises(ValueError, match="Failed attempt"):
        runtime.run_segment(plan_path, "window", 2)


def test_queue_holds_extensions_until_all_first_blocks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One worker covers every assigned seed before any second block."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "relayqueue", Path(__file__).parents[1] / "scripts/run_relay_umbrella.py"
    )
    assert spec and spec.loader
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)
    cfg = dict(
        windows=[
            dict(
                label=label,
                stages=[dict(kind="production", input=f"{label}/prod/run.mdin")],
            )
            for label in ("a", "b")
        ]
    )
    monkeypatch.setattr(
        queue,
        "verify_plan",
        lambda path: (
            dict(workers=1, coverage_gate_timeout_hours=1, production_blocks=6),
            cfg,
            tmp_path,
        ),
    )
    calls = []

    def segment(path: Path, label: str, limit: int) -> None:
        calls.append((label, limit))
        (tmp_path / label).mkdir(exist_ok=True)
        (tmp_path / label / "prod").mkdir(exist_ok=True)
        (tmp_path / label / "prod/analysis_outcome.json").write_text(
            json.dumps(dict(status="completed"))
        )
        (tmp_path / label / f"segment_receipt_limit_{limit:02d}.json").write_text(
            json.dumps(dict(production_blocks=limit, status="accepted"))
        )

    monkeypatch.setattr(queue, "run_segment", segment)
    queue.worker(tmp_path / "plan.json", 0)
    assert calls == [(label, n) for n in range(1, 7) for label in ("a", "b")]
    assert len(list((tmp_path / "claims").glob("*.json"))) == 12


def test_queue_failure_blocks_extensions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Peer failure survives as a durable failed worker record."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "relayqueuefailure", Path(__file__).parents[1] / "scripts/run_relay_umbrella.py"
    )
    assert spec and spec.loader
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)
    cfg = dict(
        windows=[
            dict(
                label=label,
                stages=[dict(kind="production", input=f"{label}/prod/run.mdin")],
            )
            for label in ("a", "b")
        ]
    )
    monkeypatch.setattr(
        queue,
        "verify_plan",
        lambda path: (
            dict(workers=2, coverage_gate_timeout_hours=1, production_blocks=6),
            cfg,
            tmp_path,
        ),
    )
    (tmp_path / "b").mkdir()
    (tmp_path / "b/window_status.json").write_text(json.dumps(dict(status="failed")))
    calls = []
    monkeypatch.setattr(
        queue, "run_segment", lambda path, label, limit: calls.append((label, limit))
    )
    with pytest.raises(RuntimeError, match="peer failed"):
        queue.worker(tmp_path / "plan.json", 0)
    assert calls == []
    assert (
        json.loads((tmp_path / "workers/worker_00/status.json").read_text())["status"]
        == "failed"
    )


def test_recovered_scc_warning_retained_late_unresolved_rejected(
    tmp_path: Path,
) -> None:
    """A later warning for an already-converged step cannot hide behind its ID."""
    warning = "QMMM SCC-DFTB: SCC-DFTB FOR STEP 2 DID NOT CONVERGE AFTER 70 cycles.\n"
    text = native_text()
    (tmp_path / "run.mdout").write_text(warning + text)
    result = native_acceptance(tmp_path, 2, [250, 350])
    assert result["electronic_history"]["recovered_trial_count"] == 1
    (tmp_path / "run.mdout").write_text(text + "\n" + warning)
    with pytest.raises(ValueError, match="electronic history"):
        native_acceptance(tmp_path, 2, [250, 350])


def test_dynamic_claim_exclusion_and_failed_claim_hold(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Concurrent workers claim different windows, then hold on explicit failure."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "relayclaims", Path(__file__).parents[1] / "scripts/run_relay_umbrella.py"
    )
    assert spec and spec.loader
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)
    cfg = dict(
        windows=[
            dict(
                label=label,
                stages=[dict(kind="production", input=f"{label}/prod/run.mdin")],
            )
            for label in ("a", "b")
        ]
    )
    monkeypatch.setattr(
        queue,
        "verify_plan",
        lambda path: (dict(workers=2, production_blocks=6), cfg, tmp_path),
    )
    path = tmp_path / "plan.json"
    assert queue.claim_task(path, 0) == ("a", 1)
    assert queue.claim_task(path, 1) == ("b", 1)
    assert queue.claim_task(path, 2) == "wait"
    (tmp_path / "a").mkdir()
    (tmp_path / "a/segment_receipt_limit_01.json").write_text(
        json.dumps(dict(status="accepted", production_blocks=1))
    )
    queue.finish_claim(path, "a", 1, None)
    # A second block cannot start while b has no first block.
    assert queue.claim_task(path, 0) == "wait"
    queue.finish_claim(path, "b", 1, "native failed")
    with pytest.raises(RuntimeError, match="Failed queue claim"):
        queue.claim_task(path, 0)


def test_analysis_failure_preserves_native_then_retries_only_derived(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A broken CSV producer leaves immutable native acceptance unchanged."""
    import sys
    from types import SimpleNamespace
    from nyl12_covalent.relay_umbrella_runtime import ensure_diagnostics

    for name in ("acceptance.json", "process.json", "restart.ncrst", "trajectory.nc"):
        (tmp_path / name).write_text("immutable-native")
    before = {p.name: digest(p) for p in tmp_path.iterdir()}

    def fail(*args: object) -> None:
        (tmp_path / "raw_features.csv").write_text("partial")
        raise ValueError("derived CSV schema failure")

    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.relay_umbrella_diagnostics",
        SimpleNamespace(block_diagnostics=fail),
    )
    assert ensure_diagnostics({}, {}, tmp_path, object()) is False
    assert (
        json.loads((tmp_path / "analysis_outcome.json").read_text())["status"]
        == "failed"
    )

    def succeed(*args: object) -> None:
        for name in (
            "raw_features.csv",
            "analysis_request.json",
            "samples_and_cross_bias.npz",
        ):
            (tmp_path / name).write_text("derived")

    monkeypatch.setitem(
        sys.modules,
        "nyl12_covalent.relay_umbrella_diagnostics",
        SimpleNamespace(block_diagnostics=succeed),
    )
    assert ensure_diagnostics({}, {}, tmp_path, object()) is True
    assert len(list(tmp_path.glob("failed_analysis_*"))) == 1
    assert all(digest(tmp_path / name) == value for name, value in before.items())


def test_restart_md_has_no_zero_step_electronic_evaluation(tmp_path: Path) -> None:
    """Continuation integrates steps1..N; cold-start step0 is not repeated."""
    text = native_text().replace(
        "QMMM SCC-DFTB: SCC-DFTB for step 0 converged in 20 cycles.\n", ""
    )
    (tmp_path / "run.mdout").write_text(text)
    assert native_acceptance(tmp_path, 2, [250, 350])[
        "every_integration_step_SCC_converged"
    ]


def test_final_equil_temperature_excludes_averages_and_rms(tmp_path: Path) -> None:
    text = (
        native_text()
        + "\n AVERAGES OVER 2 STEPS\nNSTEP = 2 TIME(PS) = .001 TEMP(K) = 180.0\nEPtot = -12\n R M S F L U C T U A T I O N S\nNSTEP = 2 TIME(PS) = .001 TEMP(K) = 50.0\nEPtot = 1.0"
    )
    (tmp_path / "run.mdout").write_text(text)
    assert native_acceptance(tmp_path, 2, [225, 375])["final_temperature_K"] == 300.0
