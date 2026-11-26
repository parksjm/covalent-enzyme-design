"""Tests for the bounded direct-US continuation adapter and scheduler."""

from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent import (
    direct_us_continuation as continuation,
    direct_us_continuation_scheduler as scheduler,
    direct_us_continuation_runtime as runtime_identity,
    direct_us_continuation_worker as worker,
    neb_umbrella_campaign as campaign,
)


class _Topology:
    def __init__(self, atom_count: int = 4165) -> None:
        self.atoms = [object()] * atom_count


class _Mask:
    def __init__(self, topology: _Topology, mask: str) -> None:
        self.topology = topology
        self.mask = mask

    def Selected(self) -> list[int]:
        return list(range(len(self.topology.atoms)))


class _Restart:
    def __init__(self, atom_count: int = 4165) -> None:
        self.coordinates = np.zeros((atom_count, 3), dtype=float)
        self.box = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])


def _synthetic_campaign(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, dict]:
    output = tmp_path / "umbrella_sampling"
    output.mkdir()
    augmented = output / "analysis/fixed_burnin_20261126/augmented"
    augmented.mkdir(parents=True)
    prepared = output / "prepared.json"
    prepared.write_text("synthetic immutable campaign\n")
    source_input_template = (
        "imin=0, irest=1, ntx=5, nstlim=10000, dt=0.0005000,\n"
        "ntt=3, gamma_ln=2.0, temp0=300.0, ig={seed},\n"
    )
    windows = []
    for index in range(63):
        label = f"window_{index + 1:03d}"
        folder = output / label
        source_stage = folder / "03_production_block_06_5ps"
        source_stage.mkdir(parents=True)
        source_stage.joinpath("run.mdin").write_text(
            source_input_template.format(seed=202800001 + index * 100 + 7)
        )
        source_stage.joinpath("umbrella.RST").write_text("synthetic restraint\n")
        restart = source_stage / "restart.ncrst"
        restart.write_bytes(f"restart {label}".encode())
        process = source_stage / "process.json"
        process.write_text('{"exit_status": 0}\n')
        acceptance = source_stage / "acceptance.json"
        acceptance.write_text(
            json.dumps(
                {
                    "checkpoint_sha256": campaign.digest(restart),
                    "completed_steps": 10000,
                    "finite_energy_temperature": True,
                    "SCC_error_detected": False,
                    "final_temperature_K": 300.0,
                }
            )
            + "\n"
        )
        folder.joinpath("window_status.json").write_text(
            json.dumps(
                {
                    "status": "completed",
                    "production_blocks": 6,
                    "production_ps": 30,
                    "completed": ["03_production_block_06_5ps"],
                }
            )
            + "\n"
        )
        windows.append(
            {
                "label": label,
                "seed": f"synthetic_seed_{label}",
                "stages": [
                    {
                        "label": stage_label,
                        "seed": 202800001 + index * 100 + stage_index,
                    }
                    for stage_index, stage_label in enumerate(
                        [
                            "01_heating_1ps",
                            "02_equilibration_5ps",
                            *[
                                f"03_production_block_{block:02d}_5ps"
                                for block in range(1, 7)
                            ],
                        ]
                    )
                ],
            }
        )
    cfg = {
        "output": str(output),
        "windows": windows,
        "topology": "synthetic_topology",
        "sander": "synthetic_sander",
        "active_mask": "@active",
        "temperature_K": 300,
        "temperature_gate_K": [250, 350],
        "bias_convention": "Amberrk:centralE=k*delta²;nativeouterlineartails;nohalf",
        "seed_base": 202800001,
        "execution_sha256": {},
        "qm_charge": -2,
        "qm_atoms": 98,
    }
    monkeypatch.setattr(continuation.campaign, "frozen", lambda _path: cfg)
    monkeypatch.setattr(continuation, "_source_paths", lambda _path: [prepared])
    monkeypatch.setattr(
        continuation.runtime_identity,
        "capture_external_runtime",
        lambda value: {
            str(Path(value["sander"]).resolve()): "a" * 64,
            str(Path(value["topology"]).resolve()): "b" * 64,
        },
    )
    monkeypatch.setattr(
        continuation.runtime_identity, "verify_external_runtime", lambda _hashes: None
    )
    monkeypatch.setattr(continuation.parmed, "load_file", lambda _path: _Topology())
    monkeypatch.setattr(continuation, "AmberMask", _Mask)
    monkeypatch.setattr(continuation.Rst7, "open", lambda _path: _Restart())
    monkeypatch.setattr(
        continuation.segments,
        "_checkpoint",
        lambda path, baseline, movable, box: (
            {"checkpoint_sha256": campaign.digest(Path(path))},
            np.asarray(baseline).copy(),
        ),
    )
    return prepared, output, cfg


def test_replace_unique_seed_preserves_other_input_and_requires_restart() -> None:
    source = "irest=1, ntx=5, nstlim=10000, ig=202800008, temp0=300.0,\n"
    actual = continuation.replace_unique_seed(source, 202800010)
    assert actual == source.replace("ig=202800008", "ig=202800010")
    with pytest.raises(ValueError, match="exactly one"):
        continuation.replace_unique_seed("irest=1, ntx=5, ig=1, ig=2", 3)
    with pytest.raises(ValueError, match="restart coordinates"):
        continuation.replace_unique_seed("irest=0, ntx=1, ig=1", 3)
    with pytest.raises(ValueError, match="32-bit"):
        continuation.replace_unique_seed("irest=1, ntx=5, ig=1", 0)


def test_prepare_pins_63_windows_and_materializes_only_blocks_07_to_12(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, output, cfg = _synthetic_campaign(tmp_path, monkeypatch)
    control = output / "continuation_60ps_20261126"
    plan_path = continuation.prepare_continuation(prepared, control, workers=40)
    plan = json.loads(plan_path.read_text())

    assert plan["prepared_sha256"] == campaign.digest(prepared)
    assert plan["window_count"] == 63
    assert len(plan["external_runtime_sha256"]) == 2
    assert plan["engine_binary_sha256"] == "a" * 64
    assert len(plan["stages"]) == 378
    assert len(plan["windows"]) == 63
    assert plan["windows"][0]["stages"][0] == {
        "label": "03_production_block_07_5ps",
        "kind": "production",
        "block": 7,
    }
    original_seeds = {
        stage["seed"] for row in cfg["windows"] for stage in row["stages"]
    }
    extension_seeds = {record["seed"] for record in plan["stages"]}
    assert len(extension_seeds) == 378
    assert extension_seeds.isdisjoint(original_seeds)
    first = plan["stages"][0]
    stage = Path(first["stage"])
    assert campaign.digest(stage / "run.mdin") == first["mdin_sha256"]
    assert campaign.digest(stage / "umbrella.RST") == first["restraint_sha256"]
    assert (control / "plan.sha256").read_text().strip() == campaign.digest(plan_path)
    verified, _ = continuation.verify_plan(plan_path)
    assert verified["prepared_sha256"] == plan["prepared_sha256"]


def test_prepare_refuses_unaccepted_original_status_before_materialization(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, output, _ = _synthetic_campaign(tmp_path, monkeypatch)
    status = output / "window_001/window_status.json"
    value = json.loads(status.read_text())
    value["status"] = "running"
    status.write_text(json.dumps(value))
    control = output / "continuation_60ps_20261126"

    with pytest.raises(ValueError, match="exactly 30 production ps"):
        continuation.prepare_continuation(prepared, control)
    assert not control.exists()
    assert not (output / "window_001/03_production_block_07_5ps").exists()


def test_verify_plan_detects_new_input_and_original_prepared_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, output, _ = _synthetic_campaign(tmp_path, monkeypatch)
    plan_path = continuation.prepare_continuation(
        prepared, output / "continuation_60ps_20261126"
    )
    first = json.loads(plan_path.read_text())["stages"][0]
    Path(first["stage"], "run.mdin").write_text("mutated\n")
    with pytest.raises(ValueError, match="run.mdin changed"):
        continuation.verify_plan(plan_path)
    Path(first["stage"], "run.mdin").write_text(
        "imin=0, irest=1, ntx=5, nstlim=10000, dt=0.0005000,\n"
        f"ntt=3, gamma_ln=2.0, temp0=300.0, ig={first['seed']},\n"
    )
    # Rebuild the plan hash only in a disposable fixture to get past the earlier
    # input tamper and establish that the immutable campaign hash is checked.
    plan = json.loads(plan_path.read_text())
    record = plan["stages"][0]
    record["mdin_sha256"] = campaign.digest(Path(record["stage"]) / "run.mdin")
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    plan_path.with_name("plan.sha256").write_text(campaign.digest(plan_path) + "\n")
    prepared.write_text("changed prepared campaign\n")
    with pytest.raises(ValueError, match="Original prepared campaign changed"):
        continuation.verify_plan(plan_path)


def test_accepted_prefix_hashes_process_file_and_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stage = tmp_path / "window_001/03_production_block_06_5ps"
    stage.mkdir(parents=True)
    process = stage / "process.json"
    process.write_text('{"exit_status": 0}\n')
    acceptance = stage / "acceptance.json"
    acceptance.write_text('{"completed_steps": 10000}\n')
    restart = stage / "restart.ncrst"
    restart.write_bytes(b"checkpoint")
    origin = {
        "stage": str(stage),
        "restart": str(restart),
        "process_sha256": campaign.digest(process),
        "acceptance_sha256": campaign.digest(acceptance),
        "restart_sha256": campaign.digest(restart),
        "atom_count": 4165,
    }
    cfg = {"active_mask": "@active"}
    row = {"label": "window_001", "seed": "synthetic seed"}
    plan = {"prepared_input_identity": {"window_001": origin}, "stages": []}
    monkeypatch.setattr(worker, "AmberMask", _Mask)
    monkeypatch.setattr(worker.Rst7, "open", lambda _path: _Restart())
    monkeypatch.setattr(
        worker.segments, "_checkpoint", lambda path, *_args: ({}, np.zeros((4165, 3)))
    )
    assert worker._validate_accepted_prefix(plan, cfg, row, _Topology(), 7) == restart


def test_accepted_prefix_validates_continuation_block_before_using_restart(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_stage = tmp_path / "window_001/03_production_block_06_5ps"
    old_stage.mkdir(parents=True)
    old_process = old_stage / "process.json"
    old_process.write_text('{"exit_status": 0}\n')
    old_acceptance = old_stage / "acceptance.json"
    old_acceptance.write_text("{}\n")
    old_restart = old_stage / "restart.ncrst"
    old_restart.write_bytes(b"old checkpoint")
    new_stage = tmp_path / "window_001/03_production_block_07_5ps"
    new_stage.mkdir()
    process = new_stage / "process.json"
    process.write_text('{"exit_status": 0}\n')
    restart = new_stage / "restart.ncrst"
    restart.write_bytes(b"new checkpoint")
    acceptance = new_stage / "acceptance.json"
    acceptance.write_text(
        json.dumps(
            {
                "checkpoint_sha256": campaign.digest(restart),
                "completed_steps": 10000,
                "finite_energy_temperature": True,
            }
        )
        + "\n"
    )
    (new_stage / "continuation_status.json").write_text('{"status":"accepted"}\n')
    origin = {
        "stage": str(old_stage),
        "restart": str(old_restart),
        "process_sha256": campaign.digest(old_process),
        "acceptance_sha256": campaign.digest(old_acceptance),
        "restart_sha256": campaign.digest(old_restart),
        "atom_count": 4165,
    }
    plan = {
        "prepared_input_identity": {"window_001": origin},
        "stages": [{"window": "window_001", "block": 7, "stage": str(new_stage)}],
    }
    monkeypatch.setattr(worker, "AmberMask", _Mask)
    monkeypatch.setattr(worker.Rst7, "open", lambda _path: _Restart())
    monkeypatch.setattr(
        worker.segments,
        "_checkpoint",
        lambda path, *_args: (
            {"checkpoint_sha256": campaign.digest(Path(path))},
            np.zeros((4165, 3)),
        ),
    )

    assert (
        worker._validate_accepted_prefix(
            plan,
            {"active_mask": "@active"},
            {"label": "window_001", "seed": "seed"},
            _Topology(),
            8,
        )
        == restart
    )


def test_marker_state_stops_at_partial_attempt_and_respects_blocked_receipt(
    tmp_path: Path,
) -> None:
    control = tmp_path / "control"
    records = []
    for block in range(7, 13):
        stage = tmp_path / "window_001" / f"block_{block}"
        stage.mkdir(parents=True)
        (stage / "run.mdin").write_text("input")
        (stage / "umbrella.RST").write_text("restraint")
        records.append({"window": "window_001", "block": block, "stage": str(stage)})
    plan = {"control_dir": str(control), "stages": records}
    assert scheduler._window_marker_state(plan, "window_001") == (7, None)
    (Path(records[0]["stage"]) / "continuation_status.json").write_text(
        '{"status":"accepted"}'
    )
    assert scheduler._window_marker_state(plan, "window_001") == (8, None)
    (Path(records[1]["stage"]) / "run.mdout").write_text("partial")
    block, reason = scheduler._window_marker_state(plan, "window_001")
    assert block is None and "artifacts without acceptance" in reason
    receipt = control / "receipts/window_001.json"
    receipt.parent.mkdir(parents=True)
    receipt.write_text('{"reason":"manual validation stop"}')
    assert scheduler._window_marker_state(plan, "window_001") == (
        None,
        "manual validation stop",
    )


def test_window_lock_is_nonblocking_and_exclusive(tmp_path: Path) -> None:
    with scheduler._window_lock(tmp_path, "window_001") as first:
        assert first is True
        with scheduler._window_lock(tmp_path, "window_001") as second:
            assert second is False


def test_dynamic_scheduler_advances_each_window_at_block_boundaries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan_path = tmp_path / "plan.json"
    plan_path.write_text("synthetic plan\n")
    rows = [{"label": "window_001"}, {"label": "window_002"}]
    progress = {row["label"]: 7 for row in rows}
    calls = []
    plan = {"workers": 1}
    cfg = {"topology": "synthetic topology", "windows": rows}

    monkeypatch.setattr(
        scheduler.continuation, "verify_plan", lambda _path: (plan, cfg)
    )
    monkeypatch.setattr(scheduler.parmed, "load_file", lambda _path: _Topology())

    @contextmanager
    def unlocked(_control: Path, _label: str):
        yield True

    def marker_state(_plan: dict, label: str) -> tuple[int | None, str | None]:
        block = progress[label]
        return (block, None) if block <= 12 else (None, None)

    def run_stage(
        _path: Path,
        _plan: dict,
        _cfg: dict,
        row: dict,
        _top: object,
        block: int,
        _worker: int,
    ) -> None:
        calls.append((row["label"], block))
        progress[row["label"]] += 1

    monkeypatch.setattr(scheduler, "_window_lock", unlocked)
    monkeypatch.setattr(scheduler, "_window_marker_state", marker_state)
    monkeypatch.setattr(scheduler, "_run_one_stage", run_stage)
    result = scheduler.run_worker(plan_path, 0, idle_seconds=0.001)

    assert result["status"] == "completed"
    assert calls[:4] == [
        ("window_001", 7),
        ("window_002", 7),
        ("window_001", 8),
        ("window_002", 8),
    ]
    assert len(calls) == 12


def test_scheduler_claim_is_exclusive_for_each_worker_id(tmp_path: Path) -> None:
    scheduler._claim_worker(tmp_path, 3, "a" * 64)
    with pytest.raises(FileExistsError):
        scheduler._claim_worker(tmp_path, 3, "a" * 64)


def test_scheduler_rechecks_running_marker_under_lock_before_blocking(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan_path = tmp_path / "plan.json"
    plan_path.write_text("synthetic plan\n")
    plan = {"workers": 1}
    cfg = {"topology": "synthetic topology", "windows": [{"label": "window_001"}]}
    calls = 0

    monkeypatch.setattr(
        scheduler.continuation, "verify_plan", lambda _path: (plan, cfg)
    )
    monkeypatch.setattr(scheduler.parmed, "load_file", lambda _path: _Topology())

    @contextmanager
    def unlocked(_control: Path, _label: str):
        yield True

    def just_finished(_plan: dict, _label: str) -> tuple[int | None, str | None]:
        nonlocal calls
        calls += 1
        if calls == 1:
            return None, "block 7 has nonaccepted attempt state running"
        return None, None

    monkeypatch.setattr(scheduler, "_window_lock", unlocked)
    monkeypatch.setattr(scheduler, "_window_marker_state", just_finished)
    monkeypatch.setattr(scheduler.time, "sleep", lambda _seconds: None)
    result = scheduler.run_worker(plan_path, 0)

    assert result["status"] == "completed"
    assert not (tmp_path / "receipts/window_001.json").exists()


def test_external_runtime_matches_relay_pin_and_detects_parameter_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "repository"
    relay = root / runtime_identity.RELAY_PREPARED
    relay.parent.mkdir(parents=True)
    parameter_dir = tmp_path / "amber" / "venv/dat/slko/3ob-3-1"
    parameter_dir.mkdir(parents=True)
    expected = {}
    for first in runtime_identity.ELEMENTS:
        for second in runtime_identity.ELEMENTS:
            path = parameter_dir / f"{first}-{second}.skf"
            path.write_text(f"{first}-{second}\n")
            expected[str(path.resolve())] = campaign.digest(path)
    relay.write_text(json.dumps({"external_runtime_sha256": expected}))
    sander = tmp_path / "amber/venv/bin/sander"
    topology = tmp_path / "synthetic.prmtop"
    sander.parent.mkdir(parents=True)
    sander.write_text("amber engine")
    topology.write_text("topology")
    monkeypatch.setattr(runtime_identity.campaign, "ROOT", root)

    hashes = runtime_identity.capture_external_runtime(
        {"sander": str(sander), "topology": str(topology)}
    )
    assert len(hashes) == 18
    runtime_identity.verify_external_runtime(hashes)
    (parameter_dir / "O-O.skf").write_text("changed parameter\n")
    with pytest.raises(ValueError, match="Pinned external runtime file changed"):
        runtime_identity.verify_external_runtime(hashes)


@pytest.mark.parametrize("mode", ["native_failure", "postprocess_failure", "success"])
def test_worker_preserves_native_vs_derived_block_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    control = tmp_path / "control"
    (control / "workers").mkdir(parents=True)
    plan_path = control / "plan.json"
    plan_path.write_text("synthetic pinned plan\n")
    stage = tmp_path / "window_001/03_production_block_07_5ps"
    stage.mkdir(parents=True)
    (stage / "run.mdin").write_text("input")
    (stage / "umbrella.RST").write_text("restraint")
    previous = tmp_path / "window_001/03_production_block_06_5ps/restart.ncrst"
    previous.parent.mkdir(parents=True)
    previous.write_bytes(b"old checkpoint")
    record = {
        "window": "window_001",
        "block": 7,
        "label": "03_production_block_07_5ps",
        "stage": str(stage),
        "mdin_sha256": campaign.digest(stage / "run.mdin"),
        "restraint_sha256": campaign.digest(stage / "umbrella.RST"),
    }
    plan = {
        "prepared_sha256": "c" * 64,
        "source_sha256": {"synthetic source": "d" * 64},
        "prepared_input_identity": {"window_001": {"restart_sha256": "e" * 64}},
        "stages": [record],
        "augmentation_root": str(tmp_path / "augmented"),
        "scc_step_verification": "not_claimed",
    }
    cfg = {
        "output": str(tmp_path),
        "sander": "/tmp/amber/bin/sander",
        "topology": "/tmp/system.prmtop",
        "temperature_gate_K": [250, 350],
        "roles_1based": {},
    }
    row = {"label": "window_001"}
    monkeypatch.setattr(worker.continuation, "verify_plan", lambda _path: (plan, cfg))
    monkeypatch.setattr(worker, "_validate_accepted_prefix", lambda *_args: previous)

    def run_native(*_args: object, **_kwargs: object) -> Path:
        if mode == "native_failure":
            raise RuntimeError("native engine failure")
        (stage / "process.json").write_text('{"exit_status": 0}\n')
        restart = stage / "restart.ncrst"
        restart.write_bytes(b"new checkpoint")
        return restart

    monkeypatch.setattr(worker.runtime, "_run_stage", run_native)
    monkeypatch.setattr(
        worker.runtime,
        "stage_acceptance",
        lambda *_args: {
            "completed_steps": 10000,
            "finite_energy_temperature": True,
            "SCC_error_detected": False,
        },
    )
    monkeypatch.setattr(
        worker,
        "_window_baseline",
        lambda *_args: (np.zeros((4, 3)), np.ones(4, dtype=bool), np.ones(6)),
    )
    monkeypatch.setattr(
        worker.segments,
        "_checkpoint",
        lambda restart, *_args: (
            {"checkpoint_sha256": campaign.digest(restart)},
            np.zeros((4, 3)),
        ),
    )
    monkeypatch.setattr(
        worker.runtime, "observables", lambda *_args: {"cv_values": [0.0, 0.0]}
    )

    def diagnostics_stub(*args: object) -> None:
        if mode == "postprocess_failure":
            raise RuntimeError("derived diagnostic failure")
        (Path(args[2]) / "samples_and_cross_bias.npz").write_bytes(b"bias")

    monkeypatch.setattr(worker.diagnostics, "block_diagnostics", diagnostics_stub)
    monkeypatch.setattr(
        worker.direct_us_coordinates,
        "augment_stage",
        lambda *_args: {"csv_sha256": "f" * 64},
    )

    if mode == "native_failure":
        with pytest.raises(RuntimeError, match="native engine failure"):
            worker._run_one_stage(plan_path, plan, cfg, row, _Topology(4), 7, 0)
    else:
        worker._run_one_stage(plan_path, plan, cfg, row, _Topology(4), 7, 0)

    status = json.loads((stage / "continuation_status.json").read_text())
    if mode == "native_failure":
        assert status["status"] == "failed"
        assert not (stage / "acceptance.json").exists()
        return
    postprocess = json.loads((stage / "postprocess.json").read_text())
    assert status["status"] == "accepted"
    assert status["production_ps_per_window"] == 35
    assert status["continuation_ps_per_window"] == 5
    assert postprocess["status"] == (
        "failed" if mode == "postprocess_failure" else "completed"
    )
    assert postprocess["native_stage_accepted"] is True
    assert (stage / "acceptance.json").is_file()
