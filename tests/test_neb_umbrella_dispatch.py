"""Bounded static ownership and exact existing-job adoption checks."""

from pathlib import Path

import pytest

from nyl12_covalent import neb_umbrella_dispatch as d


def job(
    cfg: dict, prepared: Path, label: str, ended: bool = True, code: int = 0
) -> None:
    """Write representative original launch_job metadata and scientific receipt."""
    folder = prepared.parent / label
    (folder / "job").mkdir(parents=True, exist_ok=True)
    record = {
        "command": d.command(cfg, prepared, label),
        "input_sha256": {
            str(prepared): d.digest(prepared),
            **cfg["prepared_input_sha256"],
        },
        "ended_utc": "now" if ended else None,
        "exit_status": code,
    }
    d.write(folder / "job" / "command.json", record)
    d.write(folder / "job" / "status.json", record)
    d.write(
        folder / "window_status.json",
        {
            "status": "completed",
            "production_ps": 30,
            "completed": ["heat", "production"],
        },
    )


@pytest.fixture
def campaign(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple:
    """Create immutable host identities without scientific computations."""
    monkeypatch.setattr(d.socket, "gethostname", lambda: "violet")
    labels = [f"window_{i:03d}" for i in range(1, 64)]
    payload = tmp_path / "input"
    payload.write_text("immutable")
    cfg = {
        "priority_schedule": labels,
        "windows": [
            {"label": label, "stages": [{"label": "heat"}, {"label": "production"}]}
            for label in labels
        ],
        "python": "/usr/bin/python3",
        "runtime_script": "/old/run.py",
        "sander": "/amber/bin/sander",
        "session_prefix": "test",
        "prepared_input_sha256": {str(payload): d.digest(payload)},
        "frozen_source_sha256": {},
        "execution_sha256": {},
    }
    prepared = tmp_path / "prepared.json"
    d.write(prepared, cfg)
    (tmp_path / "prepared.sha256").write_text(d.digest(prepared))
    for index, label in enumerate(labels[:8], 1):
        folder = tmp_path / "workers" / f"worker_{index:02d}"
        folder.mkdir(parents=True)
        d.write(folder / "status.json", {"current": label})
        job(cfg, prepared, label)
    hosts = {
        "violet": {"prepared": str(prepared), "sha256": d.digest(prepared)},
        "azalea": {"prepared": str(prepared), "sha256": d.digest(prepared)},
    }
    return cfg, prepared, hosts


def make_plan(campaign: tuple) -> Path:
    """Prepare fixture plan with real retained controller snapshots."""
    _, prepared, hosts = campaign
    return d.prepare(prepared, hosts, "violet", "azalea")


def test_partition(campaign: tuple) -> None:
    cfg, _, _ = campaign
    groups = d.partition(cfg["priority_schedule"], cfg["priority_schedule"][:8])
    assert len(groups["local"]) == 8
    assert [len(tasks) for tasks in groups["remote"]].count(2) == 23
    with pytest.raises(ValueError, match="unique"):
        d.partition(cfg["priority_schedule"], ["window_001"] * 8)


def test_prepare_freeze_and_duplicate(campaign: tuple) -> None:
    path = make_plan(campaign)
    plan, _, _ = d.frozen(path, "violet")
    assert len(plan["source_sha256"]) == 4
    with pytest.raises(FileExistsError):
        make_plan(campaign)
    with pytest.raises(ValueError, match="Execution host"):
        d.frozen(path, "azalea")
    path.write_text(path.read_text() + " ")
    with pytest.raises(ValueError, match="plan changed"):
        d.frozen(path, "violet")


@pytest.mark.parametrize("change", ["duplicate", "capacity", "mode"])
def test_plan_rejects(campaign: tuple, change: str) -> None:
    plan = d.read(make_plan(campaign))
    if change == "duplicate":
        plan["assignments"]["azalea"][0][0]["label"] = "window_001"
    elif change == "capacity":
        plan["assignments"]["azalea"].pop()
    else:
        plan["assignments"]["violet"][0][0]["mode"] = "launch"
    with pytest.raises(ValueError):
        d.validate_plan(plan)


def test_prepared_identity_guards(campaign: tuple) -> None:
    cfg, prepared, hosts = campaign
    hosts["violet"]["sha256"] = "changed"
    with pytest.raises(ValueError, match="Local host identity"):
        make_plan(campaign)
    hosts["violet"]["sha256"] = d.digest(prepared)
    path = make_plan(campaign)
    Path(next(iter(cfg["prepared_input_sha256"]))).write_text("changed")
    with pytest.raises(ValueError, match="Frozen source/input"):
        d.frozen(path, "violet")
    prepared.write_text(prepared.read_text() + " ")
    with pytest.raises(ValueError, match="Host prepared identity"):
        d.frozen(path, "violet")


def test_adopt_receipts(campaign: tuple) -> None:
    cfg, prepared, _ = campaign
    label = "window_001"
    job(cfg, prepared, label, ended=False)
    assert not d.complete(cfg, prepared, label)
    job(cfg, prepared, label, code=2)
    with pytest.raises(RuntimeError, match="process failed"):
        d.complete(cfg, prepared, label)
    job(cfg, prepared, label)
    assert d.complete(cfg, prepared, label)
    receipt = prepared.parent / label / "window_status.json"
    d.write(receipt, {"status": "completed", "production_ps": 5, "completed": []})
    with pytest.raises(ValueError, match="scientific window"):
        d.complete(cfg, prepared, label)
    record = prepared.parent / label / "job" / "status.json"
    value = d.read(record)
    value["command"][-1] = "window_063"
    d.write(record, value)
    with pytest.raises(ValueError, match="identity mismatch"):
        d.verify_job(cfg, prepared, label)


def test_worker_adopts_then_launches(
    campaign: tuple, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg, prepared, _ = campaign
    path = make_plan(campaign)
    launched = []

    def fake_launch(*args: object, **kwargs: object) -> None:
        label = args[3][-1]
        launched.append(label)
        job(cfg, prepared, label)

    monkeypatch.setattr(d, "launch_job", fake_launch)
    d.worker(path, "violet", 0)
    assert launched == ["window_009"]
    status = d.read(path.parent / "violet" / "worker01" / "status.json")
    assert status["completed"] == ["window_001", "window_009"]
    with pytest.raises(FileExistsError):
        d.worker(path, "violet", 0)
    with pytest.raises(ValueError, match="index"):
        d.worker(path, "violet", -1)


@pytest.mark.parametrize("failure", ["launch", "existing"])
def test_worker_failure(
    campaign: tuple, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    _, prepared, _ = campaign
    path = make_plan(campaign)

    def fail(*args: object, **kwargs: object) -> None:
        raise RuntimeError("tmux failed")

    monkeypatch.setattr(d, "launch_job", fail)
    if failure == "existing":
        folder = prepared.parent / "window_009"
        folder.mkdir()
        d.write(folder / "window_status.json", {})
    with pytest.raises((RuntimeError, FileExistsError)):
        d.worker(path, "violet", 0)
    status = d.read(path.parent / "violet" / "worker01" / "status.json")
    assert status["status"] == "failed"
    assert status["completed"] == ["window_001"]


def test_launch_bounded(campaign: tuple, monkeypatch: pytest.MonkeyPatch) -> None:
    path = make_plan(campaign)
    calls = []
    monkeypatch.setattr(d, "launch_job", lambda *a, **kw: calls.append((a, kw)))
    d.launch(path, "violet")
    assert len(calls) == 8
    assert all(call[0][3][1] == d.read(path)["script"] for call in calls)
    monkeypatch.setattr(d.socket, "gethostname", lambda: "azalea")
    calls.clear()
    d.launch(path, "azalea")
    assert len(calls) == 24


def test_partial_start_status(campaign: tuple) -> None:
    cfg, prepared, _ = campaign
    (prepared.parent / "window_001" / "job" / "status.json").write_text("")
    assert not d.complete(cfg, prepared, "window_001")
    with pytest.raises(ValueError, match="write in progress"):
        make_plan(campaign)


def test_worker_polls(campaign: tuple, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg, prepared, _ = campaign
    path = make_plan(campaign)
    job(cfg, prepared, "window_001", ended=False)
    monkeypatch.setattr(d.time, "sleep", lambda _: job(cfg, prepared, "window_001"))
    monkeypatch.setattr(d, "launch_job", lambda *a, **kw: job(cfg, prepared, a[3][-1]))
    d.worker(path, "violet", 0)
    assert (
        d.read(path.parent / "violet" / "worker01" / "status.json")["status"]
        == "completed"
    )
