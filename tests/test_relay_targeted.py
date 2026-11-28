"""Bounded preparation and fake native execution for targeted relay replicas."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from nyl12_covalent import relay_targeted as prepare
from nyl12_covalent import relay_targeted_runtime as runtime
from nyl12_covalent.neb_umbrella_campaign import digest, write_json
from nyl12_covalent.neb_umbrella_inputs import active_group_cards

SOURCE = """audited
&qmmm
 qm_theory='DFTB3', qmcharge=-2, dftb_telec=0.0,
 qmmask='@1-98', scfconv=1.d-10,
/
&dlfind
 active=':1-20,25',
/
"""


@pytest.fixture
def plan(tmp_path: Path) -> Path:
    """Provide a complete local plan with explicitly pinned fixture files."""
    parent = tmp_path / "campaign"
    parent.mkdir()
    topology = parent / "system.prmtop"
    topology.write_text("accepted topology")
    source = parent / "source.mdin"
    source.write_text(SOURCE)
    engine = parent / "sander"
    engine.write_text("fake engine")
    code = parent / "code.py"
    code.write_text("frozen code")
    restart = parent / "end30.ncrst"
    restart.write_bytes(b"\x00exact binary seed\xff")
    (parent / "analysis_request.json").write_text("{}")
    prepared = parent / "prepared.json"
    write_json(
        prepared,
        dict(
            topology="system.prmtop",
            temperature_K=300,
            system=dict(active_mask=":1-20,25", qm_atom_ids_mask="@1-98"),
            physical_hamiltonian_sha256="parent Hamiltonian",
            frozen_files_sha256={"system.prmtop": digest(topology)},
        ),
    )
    prepared.with_suffix(".sha256").write_text(digest(prepared))
    plan_path = parent / "target_plan.json"
    write_json(
        plan_path,
        dict(
            parent_prepared=str(prepared),
            output=str(parent / "additional_sampling_20261126"),
            topology=str(topology),
            source_mdin=str(source),
            sander=str(engine),
            environment={"OMP_NUM_THREADS": "1"},
            runtime_sha256={str(engine): digest(engine)},
            source_manifest={str(code): digest(code)},
            targets=[
                dict(
                    label="cross021",
                    source_restart=str(restart),
                    source_restart_sha256=digest(restart),
                    feature_names=["q_acyl", "q_del"],
                    centers_A=[1.75, 1.5],
                    k=25,
                    initial_ig=123,
                    source_lineage={"window": "022", "ps": 30},
                )
            ],
        ),
    )
    return plan_path


def test_fresh_300k_changes_only_restart_control() -> None:
    args = (SOURCE, "equil", 123, ":1-20,25")
    continuous = prepare.target_stage_input(*args)
    fresh = prepare.target_stage_input(*args, fresh=True)
    assert fresh == continuous.replace("irest=1, ntx=5", "irest=0, ntx=1, tempi=300.0")
    assert "25.0" not in fresh and "TEMP0" not in fresh
    assert "nstlim=10000" in fresh and "qmshake=0" in fresh
    assert "ntc=1, ntf=1" in fresh and "ntr=0" in fresh
    assert active_group_cards(":1-20,25") in fresh
    production = prepare.target_stage_input(SOURCE, "production", 123, ":1-20,25")
    assert "irest=1, ntx=5" in production and "tempi=" not in production
    assert "qmmask='@1-98'" in production
    with pytest.raises(ValueError, match="only for equilibration"):
        prepare.target_stage_input(SOURCE, "production", 123, ":1-20,25", True)


def test_preparation_exact_seed_and_hash_guards(plan: Path) -> None:
    prepared = prepare.prepare_targets(plan)
    cfg = prepare.verify_target(prepared, digest(prepared))
    row = cfg["windows"][0]
    assert (prepared.parent / row["seed"]).read_bytes() == Path(
        row["source_restart"]
    ).read_bytes()
    assert row["source_lineage"] == {"window": "022", "ps": 30}
    assert len(row["stages"]) == 4
    assert [spec["seed"] for spec in row["stages"]] == [123, 124, 125, 126]
    assert [spec["steps"] for spec in row["stages"]] == [2000, 8000, 10000, 10000]
    assert [spec["pmf_included"] for spec in row["stages"]] == [
        False,
        False,
        True,
        True,
    ]
    for index, spec in enumerate(row["stages"]):
        text = (prepared.parent / spec["input"]).read_text()
        assert f"nstlim={spec['steps']}," in text
        assert ("tempi=300.0" in text) == (index == 0)
        assert ("irest=1, ntx=5" in text) == (index > 0)
    for spec in row["stages"]:
        assert f"ig={spec['seed']}," in (prepared.parent / spec["input"]).read_text()
    assert cfg["production_horizon_ps"] == 10
    assert cfg["automatic_primary_analysis"] is False
    assert row["bias"]["force_constants"] == [25.0, 25.0]
    assert cfg["physical_hamiltonian_sha256"] == "parent Hamiltonian"
    with pytest.raises(FileExistsError):
        prepare.prepare_targets(plan)
    with pytest.raises(ValueError, match="prepared hash"):
        prepare.verify_target(prepared, "wrong")
    frozen = prepared.parent / row["stages"][0]["input"]
    frozen.write_text("changed")
    with pytest.raises(ValueError, match="Pinned file"):
        prepare.verify_target(prepared, digest(prepared))


@pytest.mark.parametrize(
    "update,message",
    [
        ({"k": 100}, "constant"),
        ({"label": "../escape"}, "labels"),
        ({"source_restart_sha256": "bad"}, "restart hash"),
        ({"initial_ig": -1}, "Seed"),
    ],
)
def test_bad_target_rejected_before_directory(
    plan: Path, update: dict, message: str
) -> None:
    cfg = json.loads(plan.read_text())
    cfg["targets"][0].update(update)
    write_json(plan, cfg)
    with pytest.raises(ValueError, match=message):
        prepare.prepare_targets(plan)
    assert not Path(cfg["output"]).exists()


def test_empty_hash_manifest_and_relative_path_fail(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="required"):
        prepare.verify_hashes({})
    with pytest.raises(ValueError, match="relative"):
        prepare.verify_hashes({"relative": "bad"})
    with pytest.raises(ValueError, match="Absolute"):
        prepare._absolute("relative")


def _fake_runtime(monkeypatch: pytest.MonkeyPatch) -> tuple[Mock, Mock]:
    initial = SimpleNamespace(coordinates=np.zeros((104922, 3)), box=np.arange(6))
    monkeypatch.setattr(runtime.parmed, "load_file", lambda path: object())
    monkeypatch.setattr(runtime.Rst7, "open", lambda path: initial)
    monkeypatch.setattr(
        runtime,
        "AmberMask",
        lambda top, mask: SimpleNamespace(
            Selected=lambda: range(98 if mask == "@1-98" else 4165)
        ),
    )

    def run(stage: Path, engine: Path, topology: Path, start: Path, env: dict) -> Path:
        (stage / "restart.ncrst").write_bytes(start.read_bytes() + b"next")
        return stage / "restart.ncrst"

    monkeypatch.setattr(runtime, "_run_stage", run)
    monkeypatch.setattr(runtime, "native_acceptance", lambda *args: {"accepted": True})
    checkpoint = Mock(
        return_value=({"checkpoint_sha256": "checkpoint"}, initial.coordinates)
    )
    monkeypatch.setattr(runtime, "_checkpoint", checkpoint)
    diagnostics = Mock(return_value=True)
    monkeypatch.setattr(runtime, "ensure_diagnostics", diagnostics)
    return checkpoint, diagnostics


def test_fake_pipeline_preserves_seed_boundary_and_uncensored_diagnostics(
    plan: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared = prepare.prepare_targets(plan)
    checkpoint, diagnostics = _fake_runtime(monkeypatch)
    runtime.run_target(prepared, "cross021", digest(prepared))
    status = json.loads((prepared.parent / "cross021/window_status.json").read_text())
    assert status["status"] == "completed" and status["production_ps"] == 10
    assert checkpoint.call_count == 4 and diagnostics.call_count == 4
    assert [row["pmf_included"] for row in status["stages"]] == [
        False,
        False,
        True,
        True,
    ]
    roles = [
        json.loads((call.args[2] / "sampling_role.json").read_text())
        for call in diagnostics.call_args_list
    ]
    assert [role["pmf_included"] for role in roles] == [False, False, True, True]
    for call in checkpoint.call_args_list:
        assert np.array_equal(call.args[1], np.zeros((104922, 3)))
        assert call.args[2].sum() == 4165
        assert np.array_equal(call.args[3], np.arange(6))
    # Diagnostics receive every accepted production block with no owner predicate.
    assert all(
        call.args[1]["label"] == "cross021" for call in diagnostics.call_args_list
    )
    with pytest.raises(FileExistsError, match="already attempted"):
        runtime.run_target(prepared, "cross021", digest(prepared))


def test_fake_failure_retained_and_unknown_label_rejected(
    plan: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared = prepare.prepare_targets(plan)
    _fake_runtime(monkeypatch)
    with pytest.raises(ValueError, match="exactly one"):
        runtime.run_target(prepared, "absent", digest(prepared))
    monkeypatch.setattr(
        runtime, "native_acceptance", Mock(side_effect=ValueError("SCC failed"))
    )
    with pytest.raises(ValueError, match="SCC failed"):
        runtime.run_target(prepared, "cross021", digest(prepared))
    status = json.loads((prepared.parent / "cross021/window_status.json").read_text())
    assert status["status"] == "failed" and "SCC failed" in status["error"]
    assert (prepared.parent / "cross021/01_equilibration_1ps/restart.ncrst").exists()


@pytest.mark.parametrize("change", ["box", "exterior", "velocity"])
def test_real_checkpoint_rejects_seed_boundary_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change: str
) -> None:
    from nyl12_covalent import neb_umbrella_segments as segments

    path = tmp_path / "restart.ncrst"
    path.write_bytes(b"checkpoint")
    baseline = np.zeros((3, 3))
    xyz = baseline.copy()
    box = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    changed_box = box.copy()
    velocities = np.zeros((3, 3))
    if change == "box":
        changed_box[0] += 1
    elif change == "exterior":
        xyz[2, 0] += 1e-3
    else:
        velocities[0, 0] = np.nan
    restart = SimpleNamespace(
        coordinates=xyz, velocities=velocities, hasvels=True, box=changed_box
    )
    monkeypatch.setattr(segments.Rst7, "open", lambda name: restart)
    with pytest.raises(ValueError, match="box|exterior|nonfinite"):
        segments._checkpoint(path, baseline, np.array([True, True, False]), box)


@pytest.mark.parametrize("initial", [0, -1, True, 2147483645, 2147483646, 2147483647])
def test_stage_seed_range_rejected(plan: Path, initial: int) -> None:
    cfg = json.loads(plan.read_text())
    cfg["targets"][0]["initial_ig"] = initial
    write_json(plan, cfg)
    with pytest.raises(ValueError, match="int32"):
        prepare.prepare_targets(plan)
    assert not Path(cfg["output"]).exists()


@pytest.mark.parametrize("offset", [0, 1, 2, 3, -1, -2, -3])
def test_overlapping_target_stage_seeds_rejected(plan: Path, offset: int) -> None:
    cfg = json.loads(plan.read_text())
    other = dict(cfg["targets"][0], label="cross022", initial_ig=123 + offset)
    cfg["targets"].append(other)
    write_json(plan, cfg)
    with pytest.raises(ValueError, match="unique within and across"):
        prepare.prepare_targets(plan)
    assert not Path(cfg["output"]).exists()


def test_disjoint_target_stage_seeds_and_upper_bound(plan: Path) -> None:
    cfg = json.loads(plan.read_text())
    cfg["targets"].append(
        dict(cfg["targets"][0], label="cross022", initial_ig=2147483644)
    )
    write_json(plan, cfg)
    prepared = prepare.prepare_targets(plan)
    rows = json.loads(prepared.read_text())["windows"]
    seeds = [spec["seed"] for row in rows for spec in row["stages"]]
    assert len(set(seeds)) == 8
    assert seeds[-4:] == [2147483644, 2147483645, 2147483646, 2147483647]


@pytest.mark.parametrize("steps", [2000, 8000, 10000])
def test_actual_diagnostic_grid_accepts_equilibration_horizons(steps: int) -> None:
    from nyl12_covalent.relay_umbrella_diagnostics import _validate_trajectory_times

    times = np.arange(100, steps + 1, 100, dtype=np.float64) * 0.0005
    assert _validate_trajectory_times(times, 0, steps, 100) == list(
        range(100, steps + 1, 100)
    )
