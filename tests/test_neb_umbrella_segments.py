"""Focused block-boundary resume and scheduler-adoption tests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import neb_umbrella_campaign as campaign
from nyl12_covalent import neb_umbrella_runtime as runtime
from nyl12_covalent import neb_umbrella_segments as segments


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _setup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, dict, list[Path]]:
    prepared = tmp_path / "prepared.json"
    prepared.write_text("prepared")
    top_path = tmp_path / "topology.prmtop"
    top_path.write_text("topology")
    seed = tmp_path / "seed.ncrst"
    seed.write_text("seed")
    source = tmp_path / "campaign.py"
    source.write_text("source")
    stages = [
        {"label": name, "kind": kind}
        for name, kind in [
            ("01_heating_1ps", "heat"),
            ("02_equilibration_5ps", "equil"),
            *[(f"03_production_block_{i:02d}_5ps", "production") for i in range(1, 7)],
        ]
    ]
    row = {"label": "window_001", "seed": str(seed), "stages": stages}
    cfg = {
        "windows": [row],
        "topology": str(top_path),
        "active_mask": ":active",
        "sander": "/amber/bin/sander",
        "roles_1based": {"carbon": 1, "leaving_n": 2, "thr_og": 3, "proton": 4},
        "temperature_gate_K": [250, 350],
        "heating_temperature_gate_K": [25, 350],
    }
    for spec in stages:
        stage = tmp_path / row["label"] / spec["label"]
        stage.mkdir(parents=True)
        (stage / "run.mdin").write_text("native input")
        (stage / "umbrella.RST").write_text("native restraint")
    monkeypatch.setattr(campaign, "frozen", lambda path: cfg)
    monkeypatch.setattr(campaign, "source_paths", lambda: [source])
    xyz = np.zeros((4166, 3), dtype=float)
    fake_restart = SimpleNamespace(
        coordinates=xyz,
        velocities=xyz,
        hasvels=True,
        box=np.array([30, 30, 30, 90, 90, 90], dtype=float),
    )
    monkeypatch.setattr(
        runtime.parmed,
        "load_file",
        lambda path: SimpleNamespace(atoms=[object()] * 4166),
    )
    monkeypatch.setattr(runtime.Rst7, "open", lambda path: fake_restart)
    monkeypatch.setattr(
        runtime,
        "AmberMask",
        lambda *args: SimpleNamespace(Selected=lambda: range(4165)),
    )
    monkeypatch.setattr(runtime, "observables", lambda *args: {"all_retained": True})
    monkeypatch.setattr(segments, "parmed", runtime.parmed)
    monkeypatch.setattr(segments, "Rst7", runtime.Rst7)
    monkeypatch.setattr(segments, "AmberMask", runtime.AmberMask)
    monkeypatch.setattr(
        segments,
        "_checkpoint",
        lambda path, *args: (
            {
                "checkpoint_sha256": campaign.digest(path),
                "frozen_outer_max_displacement_A": 0.0,
            },
            xyz,
        ),
    )
    monkeypatch.setattr(segments.diagnostics, "block_diagnostics", lambda *args: None)
    monkeypatch.setattr(segments.diagnostics, "export_view", lambda *args: None)
    return prepared, cfg, [tmp_path / row["label"] / spec["label"] for spec in stages]


def _native_output(stage: Path, kind: str) -> None:
    steps = {"heat": 2000, "equil": 10000, "production": 10000}[kind]
    temp = 230 if kind == "heat" else 300
    (stage / "run.mdout").write_text(
        f"NSTEP = {steps} TIME(PS) = 5.0 TEMP(K) = {temp}\n EPtot = -100.0\n"
    )
    (stage / "restart.ncrst").write_text(f"restart {stage.name}")


def _accept(stage: Path, kind: str) -> None:
    _native_output(stage, kind)
    process = {"exit_status": 0, "ended_utc": "now"}
    (stage / "process.json").write_text(json.dumps(process))
    runtime.stage_acceptance(stage, kind, [25, 350] if kind == "heat" else [250, 350])
    campaign.write_json(
        stage / "acceptance.json",
        {"checkpoint_sha256": campaign.digest(stage / "restart.ncrst")},
    )


def _fake_run(
    stage: Path, sander: Path, topology: Path, start: Path, environment: dict
) -> Path:
    kind = (
        "heat"
        if stage.name.startswith("01_")
        else "equil"
        if stage.name.startswith("02_")
        else "production"
    )
    _native_output(stage, kind)
    (stage / "process.json").write_text('{"exit_status": 0, "ended_utc": "now"}\n')
    return stage / "restart.ncrst"


def test_new_segment_stops_at_absolute_block_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, _, stages = _setup(tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(
        runtime, "_run_stage", lambda *args: calls.append(args[0]) or _fake_run(*args)
    )

    receipt = segments.run_segment(prepared, "window_001", 1)

    assert receipt["status"] == "paused_at_block_boundary"
    assert receipt["production_blocks"] == 1
    assert receipt["production_ps"] == 5
    assert calls == stages[:3]
    status = json.loads((stages[0].parent / "window_status.json").read_text())
    assert status["completed"] == [stage.name for stage in stages[:3]]
    assert (stages[0].parent / "segment_receipt_limit_01.json").is_file()


def test_six_blocks_complete_and_previous_status_is_preserved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, _, _ = _setup(tmp_path, monkeypatch)
    monkeypatch.setattr(runtime, "_run_stage", _fake_run)
    folder = prepared.parent / "window_001"

    first = segments.run_segment(prepared, "window_001", 1)
    prior_status = (folder / "window_status.json").read_bytes()
    second = segments.run_segment(prepared, "window_001", 2)
    assert (
        folder / "coverage_first_previous_window_status.json"
    ).read_bytes() == prior_status
    assert second["status"] == "paused_at_block_boundary"
    for limit in range(3, 7):
        receipt = segments.run_segment(prepared, "window_001", limit)

    assert first["production_ps"] == 5
    assert receipt["status"] == "completed"
    assert receipt["production_ps"] == 30


def test_accepted_prefix_is_skipped_and_next_restart_is_used(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, _, stages = _setup(tmp_path, monkeypatch)
    for stage, kind in zip(stages[:3], ["heat", "equil", "production"], strict=True):
        _accept(stage, kind)
    calls = []

    def run(stage: Path, *args: object) -> Path:
        calls.append((stage, args[2]))
        return _fake_run(stage, *args)

    monkeypatch.setattr(runtime, "_run_stage", run)
    receipt = segments.run_segment(prepared, "window_001", 2)

    assert calls == [(stages[3], stages[2] / "restart.ncrst")]
    assert receipt["completed"] == [stage.name for stage in stages[:4]]
    assert receipt["production_blocks"] == 2


def test_scheduler_adopts_finished_stage_and_preserves_original_process(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, _, stages = _setup(tmp_path, monkeypatch)
    for stage, kind in zip(stages[:2], ["heat", "equil"], strict=True):
        _accept(stage, kind)
    stage = stages[2]
    _native_output(stage, "production")
    original = {
        "command": ["sander", "-O"],
        "started_utc": "started",
        "cwd": str(stage),
    }
    process_path = stage / "process.json"
    process_path.write_text(json.dumps(original) + "\n")
    adoption = tmp_path / "adoption.json"
    adoption.write_text(
        json.dumps(
            {
                "stage": str(stage.resolve()),
                "native_exit_status": 0,
                "native_process_identity": {
                    "argv": original["command"],
                    "cwd": str(stage),
                    "pid": 123,
                    "start_ticks": 55,
                },
                "process_json_sha256": _digest(process_path),
                "native_outputs_complete": True,
                "observed_utc": "observed",
                "raw_native_wait_status": 0,
            }
        )
    )
    monkeypatch.setattr(
        runtime, "_run_stage", lambda *args: pytest.fail("adopt must not run Sander")
    )

    receipt = segments.run_segment(prepared, "window_001", 1, adoption)

    archived = stage / "process_before_controller_handoff.json"
    assert archived.read_text() == json.dumps(original) + "\n"
    finalized = json.loads(process_path.read_text())
    assert finalized["exit_status"] == 0
    assert finalized["ended_utc"] == "observed"
    assert finalized["completion_source"] == "verified_linux_zombie_exit_status"
    assert finalized["adoption_receipt_sha256"] == _digest(adoption)
    assert receipt["status"] == "paused_at_block_boundary"


@pytest.mark.parametrize("native_exit_status", [1, 0])
def test_adoption_mismatch_or_nonzero_is_rejected_without_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, native_exit_status: int
) -> None:
    prepared, _, stages = _setup(tmp_path, monkeypatch)
    for stage, kind in zip(stages[:2], ["heat", "equil"], strict=True):
        _accept(stage, kind)
    stage = stages[2]
    process_path = stage / "process.json"
    process_path.write_text('{"started_utc":"original"}\n')
    before = process_path.read_bytes()
    adoption = tmp_path / "adoption.json"
    adoption.write_text(
        json.dumps(
            {
                "stage": str(stage.resolve())
                if native_exit_status
                else str(tmp_path / "wrong"),
                "native_exit_status": native_exit_status,
                "native_process_identity": {
                    "argv": ["sander"],
                    "cwd": str(stage),
                    "pid": 123,
                    "start_ticks": 55,
                },
                "process_json_sha256": _digest(process_path),
                "native_outputs_complete": True,
                "observed_utc": "observed",
                "raw_native_wait_status": native_exit_status,
            }
        )
    )

    with pytest.raises(ValueError):
        segments.run_segment(prepared, "window_001", 1, adoption)

    assert process_path.read_bytes() == before
    assert not (stage / "acceptance.json").exists()


def test_rejects_second_attempt_without_partial_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, _, stages = _setup(tmp_path, monkeypatch)
    stages[0].joinpath("process.json").write_text('{"exit_status": 1}\n')
    before = (stages[0] / "process.json").read_bytes()
    monkeypatch.setattr(
        runtime, "_run_stage", lambda *args: pytest.fail("must not rerun")
    )

    with pytest.raises(FileExistsError):
        segments.run_segment(prepared, "window_001", 1)

    assert (stages[0] / "process.json").read_bytes() == before
    assert not (stages[0] / "acceptance.json").exists()


@pytest.mark.parametrize("limit", [0, 7, True])
def test_production_limit_is_bounded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, limit: int
) -> None:
    prepared, _, _ = _setup(tmp_path, monkeypatch)
    with pytest.raises(ValueError):
        segments.run_segment(prepared, "window_001", limit)


def test_zero_kelvin_seed_does_not_require_velocities(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared, cfg, _ = _setup(tmp_path, monkeypatch)
    seed = SimpleNamespace(
        coordinates=np.zeros((4166, 3)),
        velocities=np.full((4166, 3), np.nan),
        hasvels=False,
        box=np.array([30, 30, 30, 90, 90, 90], dtype=float),
    )
    monkeypatch.setattr(segments.Rst7, "open", lambda path: seed)

    segments._validate_seed(cfg, cfg["windows"][0])


def test_checkpoint_requires_finite_md_velocity_and_preserves_frozen_atoms(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "restart.ncrst"
    path.write_text("restart")
    xyz = np.zeros((2, 3))
    restart = SimpleNamespace(
        coordinates=xyz,
        velocities=xyz,
        hasvels=True,
        box=np.array([30, 30, 30, 90, 90, 90], dtype=float),
    )
    monkeypatch.setattr(segments.Rst7, "open", lambda source: restart)
    record, coordinates = segments._checkpoint(
        path, xyz, np.array([True, False]), restart.box
    )
    assert record["checkpoint_sha256"] == _digest(path)
    assert record["frozen_outer_max_displacement_A"] == 0
    assert np.array_equal(coordinates, xyz)

    restart.hasvels = False
    with pytest.raises(ValueError, match="velocities"):
        segments._checkpoint(path, xyz, np.array([True, False]), restart.box)
