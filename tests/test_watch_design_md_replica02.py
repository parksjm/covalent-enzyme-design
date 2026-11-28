"""Exercise bounded completed-replica queue without launching MD or analysis."""

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import watch_design_md_replica02 as queue


@pytest.fixture
def setup(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    names = [
        "scripts/analyze_design_md_geometry.py",
        "scripts/watch_design_md_replica02.py",
        "scripts/record_agent_action.py",
        "scripts/run_nyl12_job.py",
        "src/nyl12_covalent/jobs.py",
        "src/nyl12_covalent/__init__.py",
        "src/nyl12_covalent/design_md_geometry.py",
        "src/nyl12_covalent/design_md_completion.py",
        "src/nyl12_covalent/active_site_geometry.py",
        "src/nyl12_covalent/catalytic_alignment.py",
        "src/nyl12_covalent/catalytic_alignment_geometry.py",
        "src/nyl12_covalent/native_thr_policy.py",
        "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation_source_snapshot_checked/scripts/screen_pa6_binding_orientation.py",
    ]
    sources = []
    for name in names:
        path = project / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
        sources.append(path)
    sources.append(
        Path(
            "/data/jpc/enzyme-substrate-modeling/scripts/analyze_nyl12_nma_four_site_geometry.py"
        )
    )
    folder = tmp_path / "queue"
    plan = queue.prepare(
        project,
        folder,
        Path(sys.executable),
        ["WT", "V111C"],
        sources,
        project,
        sources[-1],
    )
    return project, folder, plan


def marker(project, case, completed=True):
    path = (
        project / queue.FOLLOW / "md" / case / "replica_02/100ns_process_complete.json"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"completed": completed}))
    return path


def test_source_pinning_and_no_overwrite(setup):
    project, folder, plan = setup
    queue.verify_sources(folder, plan)
    with pytest.raises(FileExistsError):
        queue.prepare(
            project,
            folder,
            Path(sys.executable),
            ["WT", "V111C"],
            [Path(r["path"]) for r in plan["sources"]],
            project,
            Path(plan["legacy"]),
        )
    Path(plan["sources"][0]["path"]).write_text("changed")
    with pytest.raises(ValueError, match="Pinned source changed"):
        queue.verify_sources(folder, plan)


def test_readiness_pending_failed_and_ready(setup):
    project, _, _ = setup
    assert queue.readiness(project, "WT") == "pending"
    status = project / queue.FOLLOW / "md/WT/replica_02/job/status.json"
    status.parent.mkdir(parents=True)
    status.write_text(json.dumps({"ended_utc": "done", "exit_status": 1}))
    assert queue.readiness(project, "WT") == "failed"
    marker(project, "WT", False)
    with pytest.raises(ValueError, match="Invalid completed"):
        queue.readiness(project, "WT")
    marker(project, "WT", True)
    assert queue.readiness(project, "WT") == "ready"


def test_completed_case_success_and_policy_failure(setup, monkeypatch):
    project, folder, plan = setup
    marker(project, "WT")
    marker(project, "V111C")
    monkeypatch.setattr(queue, "_journal", lambda *args: None)

    def run(command, **kwargs):
        output = Path(command[-1])
        output.mkdir(parents=True)
        dataset = command[command.index("--dataset") + 1]
        summary = {
            "dataset": dataset,
            "replica": "replica_02",
            "native_thr_role": "diagnostic_only",
            "selection_policy": "Tier A/B catalytic geometry only; native and joint metrics retained as historical diagnostics",
        }
        if dataset == "V111C":
            summary["native_thr_role"] = "gate"
        (output / "summary.json").write_text(json.dumps(summary))
        (output / "sites.csv").write_text("a\n")
        (output / "sites.jsonl").write_text("{}\n")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(queue.subprocess, "run", run)
    assert queue.watch(folder, 1, 1) == {"WT": "accepted", "V111C": "failed"}
    assert (
        json.loads((folder / "WT/scientific_acceptance.json").read_text())[
            "native_thr_policy"
        ]
        == "diagnostic_only"
    )
    assert (
        "policy"
        in json.loads((folder / "V111C/scientific_acceptance.json").read_text())[
            "error"
        ]
    )


def test_no_incomplete_analysis_and_timeout(setup, monkeypatch):
    project, folder, plan = setup
    monkeypatch.setattr(
        queue, "run_case", lambda *args: pytest.fail("must not analyze pending MD")
    )
    monkeypatch.setattr(queue.time, "monotonic", iter([0.0, 3601.0]).__next__)
    assert queue.watch(folder, 1, 1) == {"WT": "timeout", "V111C": "timeout"}


def test_journal_failure_blocks_analysis(setup, monkeypatch):
    project, folder, plan = setup
    marker(project, "WT")

    def fail(*args):
        raise queue.JournalError("journal offline")

    monkeypatch.setattr(queue, "_journal", fail)
    monkeypatch.setattr(
        queue.subprocess,
        "run",
        lambda *args, **kwargs: pytest.fail("unlogged dispatch"),
    )
    with pytest.raises(queue.JournalError):
        queue.run_case(folder, plan, "WT")


def test_prepare_rejects_missing_source(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    p = project / "foo.py"
    p.write_text("x")
    with pytest.raises(ValueError, match="Missing required"):
        queue.prepare(
            project,
            tmp_path / "queue",
            Path(sys.executable),
            ["WT", "V111C"],
            [p],
            project,
            Path(
                "/data/jpc/enzyme-substrate-modeling/scripts/analyze_nyl12_nma_four_site_geometry.py"
            ),
        )


def test_journal_calls_and_failure(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        queue.subprocess, "run", lambda command, check: calls.append(command)
    )
    queue._journal(tmp_path, tmp_path, "intent", "analyze", (tmp_path / "source",))
    assert "--file" in calls[0] and "intent" in calls[0]

    def fail(*args, **kwargs):
        raise OSError("offline")

    monkeypatch.setattr(queue.subprocess, "run", fail)
    with pytest.raises(queue.JournalError):
        queue._journal(tmp_path, tmp_path, "intent", "analyze", ())


def test_invalid_queue_parameters(tmp_path):
    with pytest.raises(ValueError, match="positive"):
        queue.watch(tmp_path, 0, 1)
    with pytest.raises(ValueError, match="Absolute"):
        queue.prepare(
            Path("relative"),
            tmp_path / "queue",
            Path(sys.executable),
            ["WT", "V111C"],
            [],
            tmp_path,
            tmp_path / "legacy",
        )


def test_failed_native_case_does_not_block_other(setup, monkeypatch):
    project, folder, plan = setup
    marker(project, "WT", False)
    marker(project, "V111C", True)
    monkeypatch.setattr(
        queue, "run_case", lambda queue_path, plan_data, dataset: "accepted"
    )
    assert queue.watch(folder, 1, 1) == {"WT": "failed", "V111C": "accepted"}
    assert (folder / "WT/scientific_acceptance.json").exists()


def test_failed_process_still_verifies_source_hashes(setup, monkeypatch):
    project, folder, plan = setup
    marker(project, "WT")
    monkeypatch.setattr(queue, "_journal", lambda *args: None)

    def mutate(command, **kwargs):
        Path(plan["sources"][0]["path"]).write_text("changed during execution")
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(queue.subprocess, "run", mutate)
    assert queue.run_case(folder, plan, "WT") == "failed"
    result = json.loads((folder / "WT/scientific_acceptance.json").read_text())
    assert "Pinned source changed" in result["error"]
