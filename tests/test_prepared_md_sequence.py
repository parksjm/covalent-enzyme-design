"""Guard sequential GPU dispatch without duplicating the scientific MD driver."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import run_prepared_md_sequence as queue


def fixture_plan(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    config = tmp_path / "md.json"
    config.write_text("{}")
    source = tmp_path / "source.py"
    source.write_text("pinned")
    output = tmp_path / "md"
    output.mkdir()
    row = {"name": "replica_01", "seed": 100, "gpu": "0"}
    record = {"source_sha256": {str(source): queue.digest(source)}}
    runner = SimpleNamespace(
        __file__=str(source),
        frozen=lambda *a: (
            {"dataset": "design", "session_prefix": "test", "replicas": [row]},
            output,
            record,
        ),
        replicas=lambda *a: [row],
        require_gpu=lambda *a: None,
    )
    monkeypatch.setattr(queue, "runtime", lambda *a: runner)
    for key in (
        "AMBERHOME",
        "LD_LIBRARY_PATH",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        monkeypatch.setenv(key, "1")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    path = tmp_path / "queue.json"
    path.write_text(
        json.dumps(
            {
                "gpu": "0",
                "source_sha256": {str(source): queue.digest(source)},
                "tasks": [{"config": str(config), "replica": "replica_01"}],
            }
        )
    )
    return path, output


def test_frozen_plan_and_identity_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, _ = fixture_plan(tmp_path, monkeypatch)
    assert len(queue.checked_tasks(path, tmp_path, "0")[1]) == 1
    with pytest.raises(ValueError, match="agree"):
        queue.checked_tasks(path, tmp_path, "1")
    doc = json.loads(path.read_text())
    doc["tasks"] *= 2
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="duplicate"):
        queue.checked_tasks(path, tmp_path, "0")
    doc["tasks"] = []
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="at least"):
        queue.checked_tasks(path, tmp_path, "0")
    Path(next(iter(doc["source_sha256"]))).write_text("changed")
    with pytest.raises(ValueError, match="source changed"):
        queue.checked_tasks(path, tmp_path, "0")


@pytest.mark.parametrize("process_code,receipt_ok", [(0, True), (1, True), (0, False)])
def test_queue_waits_for_process_and_acceptance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    process_code: int,
    receipt_ok: bool,
) -> None:
    path, output = fixture_plan(tmp_path, monkeypatch)
    monkeypatch.setattr(queue, "require_idle_gpu", lambda *a: None)
    slept = []

    def launch(session: str, job: Path, *args: object, **kwargs: object) -> None:
        job.mkdir(parents=True)
        (job / "status.json").write_text(
            json.dumps({"ended_utc": None, "exit_status": None})
        )
        (job.parent / "100ns_process_complete.json").write_text(
            json.dumps({"completed": receipt_ok})
        )

    def sleep(seconds: int) -> None:
        slept.append(seconds)
        (output / "replica_01/job/status.json").write_text(
            json.dumps(
                {
                    "ended_utc": "done",
                    "exit_status": process_code,
                }
            )
        )

    monkeypatch.setattr(queue, "launch_job", launch)
    monkeypatch.setattr(queue.time, "sleep", sleep)
    status = tmp_path / "status.json"
    if process_code or not receipt_ok:
        with pytest.raises((RuntimeError, ValueError)):
            queue.run_queue(path, tmp_path, "0", status)
        assert json.loads(status.read_text())["status"] == "failed"
    else:
        queue.run_queue(path, tmp_path, "0", status)
        assert json.loads(status.read_text())["completed"] == 1
    assert slept == [60]


def test_gpu_contention_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        queue.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout="1234\n")
    )
    with pytest.raises(RuntimeError, match="existing compute"):
        queue.require_idle_gpu("0")
    monkeypatch.setattr(
        queue.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout="")
    )
    queue.require_idle_gpu("0")


def test_cli_uses_explicit_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, _ = fixture_plan(tmp_path, monkeypatch)
    seen = []
    monkeypatch.setattr(queue, "run_queue", lambda *a: seen.append(a))
    monkeypatch.setattr(
        queue.sys,
        "argv",
        [
            "run",
            "--config",
            str(path),
            "--repo",
            str(tmp_path),
            "--gpu",
            "0",
            "--status",
            str(tmp_path / "s.json"),
        ],
    )
    queue.main()
    assert seen[0][2] == "0"
