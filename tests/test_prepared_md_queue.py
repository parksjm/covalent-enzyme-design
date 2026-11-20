"""Balanced serial GPU queues and scientific-error stop receipts, no engine launch."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from nyl12_covalent import prepared_md_queue_plan as helper
from scripts import queue_prepared_md as queue

ROOT = Path(__file__).resolve().parents[1]


def fixture_queue(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, SimpleNamespace]:
    configs = []
    records = {}
    for i, name in enumerate(("neutral_asp", "neutral_ash", "charged_asp")):
        folder = tmp_path / name
        folder.mkdir()
        path = folder / "config.json"
        output = folder / "runs"
        output.mkdir()
        cfg = {
            "dataset": name,
            "session_prefix": name,
            "replicas": [
                {
                    "name": f"replica_{j + 1:02d}",
                    "seed": 100 + i * 100 + j,
                    "gpu": str(j % 2),
                }
                for j in range(4)
            ],
        }
        path.write_text(json.dumps(cfg))
        (output / "plan.json").write_text("{}")
        for row in cfg["replicas"]:
            (output / row["name"]).mkdir()
        records[str(path)] = (
            cfg,
            output,
            {"source_sha256": {str(path): helper.digest(path)}},
        )
        configs.append(str(path))
    fake = SimpleNamespace(
        __file__=str(ROOT / "student_handoff/nyl12_ash274_md/tools/run_prepared_md.py"),
        plan=lambda *a: None,
        frozen=lambda p, repo: records[str(p)],
        replicas=lambda cfg: cfg["replicas"],
        engines=lambda *a: (
            Path("/explicit/amber"),
            Path("/explicit/amber/bin/pmemd.cuda"),
        ),
        require_gpu=lambda gpu: {"gpu_index": gpu},
    )
    monkeypatch.setattr(helper, "runtime", lambda repo: fake)
    monkeypatch.setattr(queue, "runtime", lambda repo: fake)
    path = tmp_path / "queue.json"
    path.write_text(
        json.dumps(
            {
                "configs": configs,
                "gpus": ["0", "1"],
                "output": "queue",
                "session_prefix": "matched",
            }
        )
    )
    return path, fake


def test_balanced_plan_and_frozen_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, _ = fixture_queue(tmp_path, monkeypatch)
    output = helper.plan_queue(path, ROOT, existing_plans=True)
    _, record = helper.frozen_queue(path, ROOT)
    assert record["total_jobs"] == 12 and record["unique_stage_seeds"] == 72
    assert [t["dataset"] for t in record["tasks"]["0"]][:3] == [
        "neutral_asp",
        "neutral_ash",
        "charged_asp",
    ]
    assert [t["dataset"] for t in record["tasks"]["1"]][:3] == [
        "neutral_ash",
        "charged_asp",
        "neutral_asp",
    ]
    assert all(len(record["tasks"][g]) == 6 for g in ("0", "1"))
    assert not list(output.glob("worker*"))
    first = Path(record["tasks"]["0"][0]["config"])
    first.write_text("changed")
    with pytest.raises(ValueError, match="changed"):
        helper.frozen_queue(path, ROOT)


def test_seed_overlap_and_gpu_counts_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, fake = fixture_queue(tmp_path, monkeypatch)
    configs = json.loads(path.read_text())["configs"]
    row = fake.frozen(Path(configs[1]), ROOT)[0]["replicas"][0]
    row["seed"] = 100
    with pytest.raises(ValueError, match="overlap"):
        helper.plan_queue(path, ROOT, True)
    row["seed"] = 200
    row["gpu"] = "1"
    with pytest.raises(ValueError, match="exactly two"):
        helper.plan_queue(path, ROOT, True)


@pytest.mark.parametrize("failure", [False, True])
def test_worker_sequential_and_stops_on_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: bool
) -> None:
    import nyl12_covalent.jobs as jobs

    path, _ = fixture_queue(tmp_path, monkeypatch)
    output = helper.plan_queue(path, ROOT, True)
    calls = []

    def launch(
        session: str, job: Path, cwd: Path, command: list, **kwargs: object
    ) -> None:
        calls.append(command)
        job.mkdir()
        (job / "status.json").write_text(
            json.dumps({"ended_utc": "now", "exit_status": 1 if failure else 0}) + "\n"
        )
        (job.parent / "100ns_process_complete.json").write_text('{"completed":true}')

    monkeypatch.setattr(jobs, "launch_job", launch)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    if failure:
        with pytest.raises(RuntimeError, match="Replica failed"):
            queue.worker(path, ROOT, "0")
        assert len(calls) == 1
        assert (
            json.loads((output / "worker_gpu0/status.json").read_text())["status"]
            == "failed"
        )
    else:
        queue.worker(path, ROOT, "0")
        assert len(calls) == 6
        assert (
            json.loads((output / "worker_gpu0/status.json").read_text())["completed"]
            == 6
        )
        with pytest.raises(FileExistsError):
            queue.worker(path, ROOT, "0")
    assert not (output / "worker_gpu1").exists()


def test_launch_only_two_workers_after_gpu_checks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import nyl12_covalent.jobs as jobs

    path, _ = fixture_queue(tmp_path, monkeypatch)
    helper.plan_queue(path, ROOT, True)
    calls = []
    monkeypatch.setattr(jobs, "launch_job", lambda *a, **k: calls.append((a, k)))
    queue.launch(path, ROOT)
    assert len(calls) == 2 and {
        c[1]["env_overrides"]["CUDA_VISIBLE_DEVICES"] for c in calls
    } == {"0", "1"}
    assert all("worker" in c[0][3] for c in calls)


def test_worker_waits_for_each_finished_child_before_next(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import nyl12_covalent.jobs as jobs

    path, _ = fixture_queue(tmp_path, monkeypatch)
    output = helper.plan_queue(path, ROOT, True)
    active = []
    completed = []

    def launch(
        session: str, job: Path, cwd: Path, command: list, **kwargs: object
    ) -> None:
        assert not active
        job.mkdir()
        (job / "status.json").write_bytes(b'{"ended_utc":null')
        active.append(job)

    def wait(seconds: int) -> None:
        assert seconds == 30 and len(active) == 1
        job = active.pop()
        completed.append(job)
        (job / "status.json").write_text('{"ended_utc":"now","exit_status":0}\n')
        (job.parent / "100ns_process_complete.json").write_text('{"completed":true}')

    monkeypatch.setattr(jobs, "launch_job", launch)
    monkeypatch.setattr(queue.time, "sleep", wait)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    queue.worker(path, ROOT, "0")
    assert len(completed) == 6 and not active
    assert (
        json.loads((output / "worker_gpu0/status.json").read_text())["status"]
        == "completed"
    )


def test_cli_plan_never_launches_and_worker_requires_explicit_gpu(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, _ = fixture_queue(tmp_path, monkeypatch)
    monkeypatch.setattr(
        queue.sys,
        "argv",
        [
            "queue",
            "plan",
            "--config",
            str(path),
            "--repo",
            str(ROOT),
            "--existing-plans",
        ],
    )
    queue.main()
    assert (tmp_path / "queue/queue_plan.json").is_file()
    assert not list((tmp_path / "queue").glob("worker*"))
    monkeypatch.setattr(
        queue.sys,
        "argv",
        ["queue", "worker", "--config", str(path), "--repo", str(ROOT)],
    )
    with pytest.raises(SystemExit):
        queue.main()
