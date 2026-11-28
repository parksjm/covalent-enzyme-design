"""Source-acceptance and resource-guard tests without running an MSA engine."""

import importlib.util
import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace
from typing import Any

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot"
    / "wt_paired_precursor_controls_20261127/precursor_search/run_local_search.py"
)


@pytest.fixture
def runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    spec = importlib.util.spec_from_file_location("local_precursor_runner", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "HERE", tmp_path)
    db = tmp_path / "database"
    db.mkdir()
    for name in ["uniref30_2302_db", "colabfold_envdb_202108_db"]:
        (db / name).write_text("db")
    monkeypatch.setattr(module, "DATABASE", db)
    (tmp_path / "nyl12_wt_precursor.fasta").write_text(
        ">query\n" + "A" * 232 + "C" * 96 + "\n"
    )
    memory = tmp_path / "meminfo"
    memory.write_text("MemAvailable: 268435456 kB\n")
    monkeypatch.setattr(module, "Path", lambda _: memory)
    monkeypatch.setattr(
        module.shutil, "disk_usage", lambda _: SimpleNamespace(free=10**11)
    )
    monkeypatch.setattr(module, "_journal", lambda *_: None)
    monkeypatch.setattr(module, "_gpu_idle", lambda: True)
    return module


def _native_result(
    runner: Any, monkeypatch: pytest.MonkeyPatch, code: int, *, valid: bool = True
) -> None:
    def execute(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess:
        assert kwargs["env"]["CUDA_VISIBLE_DEVICES"] == runner.GPU
        assert command[command.index("--threads") + 1] == "2"
        assert command[command.index("--use-templates") + 1] == "0"
        if code == 0:
            content = ">query\n" + "A" * 232 + "C" * 96 + "\n"
            if valid:
                content += ">same_precursor\n" + "V" * 232 + "L" * 96 + "\n"
            # Real ColabFold search names its result from the input query header.
            (runner.HERE / "local_attempt01/query.a3m").write_text(content)
        return subprocess.CompletedProcess(command, code)

    monkeypatch.setattr(runner.subprocess, "run", execute)


def test_success_waits_for_idle_then_publishes_validated_source(
    runner: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    status = iter([False, True])
    waits = []
    monkeypatch.setattr(runner, "_gpu_idle", lambda: next(status))
    monkeypatch.setattr(runner.time, "sleep", waits.append)
    _native_result(runner, monkeypatch, 0)
    runner.main()
    receipt = json.loads((runner.HERE / "source_acceptance.json").read_text())
    assert receipt["status"] == "accepted"
    assert receipt["paired_homolog_rows_after_native_dedup"] == 1
    assert receipt["cleavage_after_precursor_residue1"] == 232
    assert waits == [120]


def test_nonzero_exit_cannot_publish_scientific_acceptance(
    runner: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    _native_result(runner, monkeypatch, 7)
    with pytest.raises(RuntimeError, match="Local search failed"):
        runner.main()
    assert not (runner.HERE / "source_acceptance.json").exists()
    assert (
        json.loads((runner.HERE / "local_search_process.json").read_text())["status"]
        == "process_failed"
    )


def test_zero_exit_without_usable_homologs_is_not_accepted(
    runner: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    _native_result(runner, monkeypatch, 0, valid=False)
    with pytest.raises(ValueError, match="No independently retained"):
        runner.main()
    assert not (runner.HERE / "source_acceptance.json").exists()


@pytest.mark.parametrize("guard", ["memory", "disk"])
def test_insufficient_resources_fail_before_native_launch(
    runner: Any, monkeypatch: pytest.MonkeyPatch, guard: str
) -> None:
    if guard == "memory":
        (runner.HERE / "meminfo").write_text("MemAvailable: 1 kB\n")
    else:
        monkeypatch.setattr(
            runner.shutil, "disk_usage", lambda _: SimpleNamespace(free=1)
        )
    with pytest.raises(RuntimeError):
        runner.main()
    assert not (runner.HERE / "local_attempt01").exists()


def test_actual_gpu_query_and_journal_arguments_are_explicit(
    runner: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    spec = importlib.util.spec_from_file_location(
        "local_precursor_actual_helpers", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    calls = []

    def execute(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, stdout=f"{module.GPU}, 123\n")

    monkeypatch.setattr(module.subprocess, "run", execute)
    assert module._gpu_idle() is False
    module._journal("intent", "test")
    assert "--phase" in calls[-1] and "intent" in calls[-1]


def test_multiple_query_records_rejected_without_native_execution(runner: Any) -> None:
    (runner.HERE / "nyl12_wt_precursor.fasta").write_text(
        ">a\n" + "A" * 328 + "\n>b\n" + "A" * 328 + "\n"
    )
    with pytest.raises(ValueError, match="one328-residue"):
        runner.main()


def _recovery(tmp_path: Path) -> tuple[Any, dict]:
    spec = importlib.util.spec_from_file_location(
        "completed_precursor_acceptance", SCRIPT.with_name("accept_completed_search.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    query = tmp_path / "nyl12_wt_precursor.fasta"
    query.write_text(">named_query\n" + "A" * 232 + "C" * 96 + "\n")
    output = tmp_path / "local_attempt01"
    output.mkdir()
    (output / "named_query.a3m").write_text(
        query.read_text() + ">origin\n" + "V" * 232 + "L" * 96 + "\n"
    )
    process = {
        "exit_status": 0,
        "status": "process_complete",
        "query_sha256": hashlib.sha256(query.read_bytes()).hexdigest(),
        "command": ["colabfold_search", str(query), "database", str(output)],
    }
    (tmp_path / "local_search_process.json").write_text(json.dumps(process))
    return module, process


def test_completed_search_accepts_native_header_named_file_without_rerun(
    tmp_path: Path,
) -> None:
    module, _ = _recovery(tmp_path)
    result = module.accept_source(tmp_path)
    assert result["status"] == "accepted"
    assert Path(result["precursor_a3m"]).name == "named_query.a3m"
    assert result["paired_homolog_rows_after_native_dedup"] == 1
    with pytest.raises(FileExistsError):
        module.accept_source(tmp_path)


@pytest.mark.parametrize("failure", ["exit", "changed_query", "length", "directory"])
def test_completed_search_rejects_unproven_or_changed_sources(
    tmp_path: Path, failure: str
) -> None:
    module, process = _recovery(tmp_path)
    if failure == "exit":
        process["exit_status"] = 7
    elif failure == "changed_query":
        (tmp_path / "nyl12_wt_precursor.fasta").write_text(">changed\nAAAA\n")
    elif failure == "length":
        query = tmp_path / "nyl12_wt_precursor.fasta"
        query.write_text(">short\nAAAA\n")
        process["query_sha256"] = hashlib.sha256(query.read_bytes()).hexdigest()
    else:
        process["command"][-1] = str(tmp_path / "unrelated")
    (tmp_path / "local_search_process.json").write_text(json.dumps(process))
    with pytest.raises(ValueError):
        module.accept_source(tmp_path)
    assert not (tmp_path / "source_acceptance.json").exists()
