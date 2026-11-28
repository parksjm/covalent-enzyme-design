"""Coverage receipt and queue handoff tests for NEB-US sampling."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nyl12_covalent import neb_umbrella_coverage as coverage
from nyl12_covalent import neb_umbrella_coverage_worker as worker
from nyl12_covalent.neb_umbrella_campaign import digest


LABELS = [f"window_{index:03d}" for index in range(1, 64)]
ARTIFACTS = (
    "acceptance.json",
    "process.json",
    "analysis_request.json",
    "samples_and_cross_bias.npz",
    "observables/observables.csv",
    "observables/observables.json",
)


def accepted_block(
    root: Path, label: str, *, steps: int = 10000, code: int = 0
) -> Path:
    """Write the minimal complete first-block receipt/artifact set."""
    stage = root / label / "03_production_block_01_5ps"
    for name in ARTIFACTS:
        path = stage / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("artifact\n")
    (stage / "acceptance.json").write_text(
        json.dumps({"completed_steps": steps, "finite_energy_temperature": True})
    )
    (stage / "process.json").write_text(
        json.dumps({"exit_status": code, "ended_utc": "now"})
    )
    return stage


def test_coverage_receipt_missing_then_complete(tmp_path: Path) -> None:
    """Only all complete and accepted first blocks produce a receipt."""
    assert coverage.coverage_receipt(tmp_path, LABELS) is None
    for label in LABELS:
        accepted_block(tmp_path, label)
    receipt = coverage.coverage_receipt(tmp_path, LABELS)
    assert receipt is not None
    assert receipt["status"] == "all_windows_have_first_block"
    assert receipt["windows"] == 63
    assert receipt["first_block_acceptance_sha256"][LABELS[0]] == digest(
        tmp_path / LABELS[0] / "03_production_block_01_5ps/acceptance.json"
    )


@pytest.mark.parametrize(
    ("steps", "code"), [(9999, 0), (10000, 1)], ids=["incomplete", "failed-process"]
)
def test_coverage_receipt_rejects_unaccepted_block(
    tmp_path: Path, steps: int, code: int
) -> None:
    """A present artifact set cannot mask an incomplete or failed native block."""
    accepted_block(tmp_path, LABELS[0], steps=steps, code=code)
    with pytest.raises(ValueError):
        coverage.coverage_receipt(tmp_path, [LABELS[0]])


def _release_context(
    monkeypatch: pytest.MonkeyPatch, schedule: Path, root: Path
) -> None:
    """Supply a deterministic frozen 63-window scheduling plan to release()."""
    plan = {
        "assignments": {"violet": [[{"label": label}] for label in LABELS]},
        "windows": LABELS,
    }
    monkeypatch.setattr(
        coverage,
        "context",
        lambda *_: (plan, {}, root / "prepared.json"),
    )


def test_release_rejects_count_and_hash_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Release requires all 63 exact acceptance hashes and a valid inventory."""
    schedule = tmp_path / "schedule.json"
    schedule.write_text("schedule")
    root = tmp_path / "campaign"
    (tmp_path / "violet").mkdir()
    for label in LABELS:
        accepted_block(root, label)
    _release_context(monkeypatch, schedule, root)
    receipt = coverage.coverage_receipt(root, LABELS)
    assert receipt is not None
    with pytest.raises(ValueError, match="coverage"):
        coverage.release(schedule, "violet", {**receipt, "windows": 62})
    with pytest.raises(ValueError, match="coverage"):
        coverage.release(
            schedule, "violet", {**receipt, "first_block_acceptance_sha256": {}}
        )
    altered = dict(receipt["first_block_acceptance_sha256"])
    altered[LABELS[0]] = "0" * 64
    with pytest.raises(ValueError, match="identity"):
        coverage.release(
            schedule,
            "violet",
            {**receipt, "first_block_acceptance_sha256": altered},
        )
    coverage.release(schedule, "violet", receipt)
    assert (tmp_path / "violet/coverage_gate.json").is_file()


def test_context_validates_schedule_and_host(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Context binds the scheduling file to its digest and explicit host."""
    prepared = tmp_path / "prepared.json"
    prepared.write_text("prepared")
    schedule = tmp_path / "plan.json"
    schedule.write_text(
        json.dumps(
            {
                "source_sha256": {},
                "hosts": {
                    "violet": {"prepared": str(prepared), "sha256": digest(prepared)}
                },
            }
        )
    )
    schedule.with_suffix(".sha256").write_text(digest(schedule))
    monkeypatch.setattr(coverage.socket, "gethostname", lambda: "violet.local")
    monkeypatch.setattr(coverage, "frozen", lambda _path: {"frozen": True})
    plan, cfg, actual = coverage.context(schedule, "violet")
    assert plan["hosts"]["violet"]["prepared"] == str(prepared)
    assert cfg == {"frozen": True}
    assert actual == prepared
    with pytest.raises(ValueError, match="execution host"):
        coverage.context(schedule, "other")
    schedule.write_text(schedule.read_text() + " ")
    with pytest.raises(ValueError, match="plan changed"):
        coverage.context(schedule, "violet")


def test_event_appends_durable_action(tmp_path: Path) -> None:
    """Controller events are append-only JSON lines with supplied details."""
    coverage.event(tmp_path, "handoff", pid=42)
    coverage.event(tmp_path, "resumed", window="window_001")
    rows = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text().splitlines()
    ]
    assert [row["action"] for row in rows] == ["handoff", "resumed"]
    assert rows[0]["pid"] == 42


def test_inventory_captures_exact_existing_owners(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Inventory records exact controller, window driver, and native identities."""
    root = tmp_path / "campaign"
    root.mkdir()
    prepared = root / "prepared.json"
    prepared.write_text("prepared")
    host_folder = tmp_path / "schedule/violet"
    host_folder.mkdir(parents=True)
    schedule = tmp_path / "schedule/plan.json"
    schedule.write_text("schedule")
    schedule.with_suffix(".sha256").write_text(digest(schedule))
    original = tmp_path / "old/plan.json"
    original.parent.mkdir()
    original.write_text(json.dumps({"script": "/old/dispatch.py"}))
    claim = original.parent / "violet/worker01/claim"
    claim.parent.mkdir(parents=True)
    claim.write_text("103\n")
    native_cwd = root / LABELS[0] / "03_production_block_04_5ps"
    native_cwd.mkdir(parents=True)
    cfg = {
        "sander": "/amber/bin/sander",
        "python": "/python",
        "runtime_script": "/old/run_neb_umbrella.py",
    }
    plan = {
        "source_sha256": {},
        "hosts": {"violet": {"prepared": str(prepared), "sha256": digest(prepared)}},
        "parent_plans": {"violet": str(original)},
        "assignments": {"violet": [[{"label": LABELS[0]}, {"label": LABELS[1]}]]},
    }
    monkeypatch.setattr(coverage, "context", lambda *_: (plan, cfg, prepared))
    proc_root = tmp_path / "proc"
    proc_entry = proc_root / "101"
    proc_entry.mkdir(parents=True)
    (proc_entry / "cmdline").write_bytes(b"/amber/bin/sander\0-O\0")
    original_iterdir = Path.iterdir

    def fake_iterdir(path: Path):
        return proc_root.iterdir() if path == Path("/proc") else original_iterdir(path)

    monkeypatch.setattr(Path, "iterdir", fake_iterdir)
    identities = {
        101: {"pid": 101, "ppid": 102, "cwd": str(native_cwd), "argv": [cfg["sander"]]},
        102: {
            "pid": 102,
            "ppid": 1,
            "cwd": str(root),
            "argv": [
                cfg["python"],
                cfg["runtime_script"],
                "window",
                "--config",
                str(prepared),
                "--window",
                LABELS[0],
            ],
        },
        103: {
            "pid": 103,
            "ppid": 1,
            "cwd": str(root),
            "argv": [
                cfg["python"],
                "/old/dispatch.py",
                "worker",
                "--plan",
                str(original),
                "--host",
                "violet",
                "--worker",
                "0",
            ],
        },
    }
    monkeypatch.setattr(coverage, "capture", lambda pid: identities[pid])
    inventory = coverage.inventory(schedule, "violet")
    record = json.loads(inventory.read_text())
    assert record["workers"][0]["current"] == LABELS[0]
    assert record["workers"][0]["queued"] == [LABELS[1]]
    assert record["workers"][0]["native"]["pid"] == 101


def test_activate_retires_only_recorded_controllers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Activation validates identities, records actions, and launches replacement."""
    schedule_dir = tmp_path / "schedule"
    host_folder = schedule_dir / "violet"
    host_folder.mkdir(parents=True)
    schedule = schedule_dir / "plan.json"
    schedule.write_text("schedule")
    prepared = tmp_path / "campaign/prepared.json"
    prepared.parent.mkdir()
    prepared.write_text("prepared")
    identity = {
        "pid": 1,
        "start_ticks": 9,
        "argv": ["python"],
        "cwd": "/tmp",
        "exe": "/python",
    }
    row = {"controller": identity, "driver": identity, "native": identity}
    (host_folder / "inventory.json").write_text(
        json.dumps({"schedule_sha256": digest(schedule), "workers": [row]})
    )
    cfg = {
        "python": "/python",
        "runtime_script": "/run.py",
        "sander": "/amber/bin/sander",
        "prepared_input_sha256": {},
    }
    plan = {"script": "/coverage.py", "hosts": {"violet": {"prepared": str(prepared)}}}
    monkeypatch.setattr(coverage, "context", lambda *_: (plan, cfg, prepared))
    monkeypatch.setattr(coverage, "capture", lambda _pid: identity)
    signaled = []
    monkeypatch.setattr(
        coverage, "signal_verified", lambda item, sig: signaled.append((item, sig))
    )
    monkeypatch.setattr(
        coverage,
        "suspend_driver",
        lambda native, driver: signaled.append((native, driver)),
    )
    launched = []
    monkeypatch.setattr(
        coverage, "launch_job", lambda *args, **kwargs: launched.append((args, kwargs))
    )
    coverage.activate(schedule, "violet")
    assert len(launched) == 1
    assert (host_folder / "activation.json").is_file()
    events = [
        json.loads(line)["action"]
        for line in (host_folder / "events.jsonl").read_text().splitlines()
    ]
    assert events == [
        "retire_queue_intent",
        "retire_queue_outcome",
        "suspend_window_driver_intent",
        "suspend_window_driver_outcome",
        "replacement_worker_launched",
    ]


def _worker_fixture(tmp_path: Path) -> tuple[Path, Path, dict]:
    """Create one worker's assigned current and queued windows."""
    schedule_dir = tmp_path / "schedule"
    host_dir = schedule_dir / "violet"
    host_dir.mkdir(parents=True)
    root = tmp_path / "campaign"
    root.mkdir()
    prepared = root / "prepared.json"
    prepared.write_text("immutable")
    schedule = schedule_dir / "plan.json"
    schedule.write_text("plan")
    row = {
        "current": LABELS[0],
        "queued": [LABELS[1]],
        "native": {
            "pid": 55,
            "cwd": str(root / LABELS[0] / "03_production_block_04_5ps"),
        },
        "driver": {"pid": 54},
    }
    active_stage = Path(row["native"]["cwd"])
    active_stage.mkdir(parents=True)
    (active_stage / "process.json").write_text(json.dumps({"command": ["sander"]}))
    for block in range(1, 4):
        previous = root / LABELS[0] / f"03_production_block_{block:02d}_5ps"
        previous.mkdir(parents=True)
        (previous / "acceptance.json").write_text("accepted")
    (host_dir / "inventory.json").write_text(json.dumps({"workers": [row]}))
    (host_dir / "worker01").mkdir()
    return schedule, prepared, row


def test_worker_adopts_current_runs_queue_then_waits_for_gate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Queued first blocks precede extensions, and both windows reach six blocks."""
    schedule, prepared, row = _worker_fixture(tmp_path)
    root = prepared.parent
    monkeypatch.setattr(worker, "context", lambda *_: ({}, {}, prepared))
    monkeypatch.setattr(worker, "native_completion", lambda _identity: 0)
    monkeypatch.setattr(worker, "release_driver", lambda _identity: None)
    calls: list[tuple[str, int, bool]] = []
    gate_ready = [False]

    def run_segment(
        _prepared: Path, label: str, limit: int, adoption_path: Path | None = None
    ) -> None:
        if limit > 1 and adoption_path is None:
            assert gate_ready[0]
        calls.append((label, limit, adoption_path is not None))
        stage = root / label / f"03_production_block_{limit:02d}_5ps"
        stage.mkdir(parents=True, exist_ok=True)
        (stage / "acceptance.json").write_text("accepted")

    monkeypatch.setattr(worker, "run_segment", run_segment)
    slept = []

    def create_global_gate(_seconds: int | float) -> None:
        slept.append(True)
        gate_ready[0] = True
        gate = schedule.parent / "violet/coverage_gate.json"
        gate.write_text(
            json.dumps({"status": "all_windows_have_first_block", "windows": 63})
        )

    monkeypatch.setattr(worker.time, "sleep", create_global_gate)
    worker.worker(schedule, "violet", 0, interval=1)
    assert calls[:2] == [(LABELS[0], 4, True), (LABELS[1], 1, False)]
    assert slept
    assert gate_ready[0]
    assert all(
        any(label == assigned and limit == 6 for label, limit, _ in calls)
        for assigned in (LABELS[0], LABELS[1])
    )


def test_worker_native_failure_does_not_start_queued_window(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed original native block stops the worker before segment execution."""
    schedule, prepared, _ = _worker_fixture(tmp_path)
    monkeypatch.setattr(worker, "context", lambda *_: ({}, {}, prepared))
    monkeypatch.setattr(worker, "native_completion", lambda _identity: 7)
    called = []
    monkeypatch.setattr(
        worker, "run_segment", lambda *args, **kwargs: called.append(args)
    )
    with pytest.raises(RuntimeError, match="Original native block failed"):
        worker.worker(schedule, "violet", 0, interval=1)
    assert called == []


def test_worker_rejects_invalid_global_gate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A malformed release cannot authorize extension sampling."""
    schedule, prepared, _ = _worker_fixture(tmp_path)
    root = prepared.parent
    monkeypatch.setattr(worker, "context", lambda *_: ({}, {}, prepared))
    monkeypatch.setattr(worker, "native_completion", lambda _identity: 0)
    monkeypatch.setattr(worker, "release_driver", lambda _identity: None)
    calls = []

    def run_segment(
        _prepared: Path, label: str, limit: int, adoption_path: Path | None = None
    ) -> None:
        calls.append(limit)
        stage = root / label / f"03_production_block_{limit:02d}_5ps"
        stage.mkdir(parents=True, exist_ok=True)
        (stage / "acceptance.json").write_text("accepted")

    monkeypatch.setattr(worker, "run_segment", run_segment)
    (schedule.parent / "violet/coverage_gate.json").write_text(
        json.dumps({"status": "all_windows_have_first_block", "windows": 62})
    )
    with pytest.raises(ValueError, match="Invalid coverage release"):
        worker.worker(schedule, "violet", 0, interval=1)
    assert calls == [4, 1]
