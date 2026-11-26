"""Exercise flat native producer receipts and conditional hydration gates."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest

from nyl12_covalent import hydration_followon_gates as gates
from nyl12_covalent import hydration_followon_queue as queue


def put(path: Path, value: Any) -> None:
    """Write a JSON fixture."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def native(tmp_path: Path) -> Path:
    """Use actual producer-shaped native flat acceptance and process records."""
    original = Path(
        "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
    )
    block = tmp_path / "extend_01"
    block.mkdir()
    for name in ("status.json", "scientific_acceptance.json", "process.json"):
        put(block / name, gates.read(original / "extend_01" / name))
    previous = tmp_path / "equil/restart.ncrst"
    previous.parent.mkdir()
    previous.write_bytes(b"predecessor")
    inputs = block / "inputs"
    inputs.mkdir()
    (inputs / "restart.ncrst").write_bytes(previous.read_bytes())
    (block / "restart.ncrst").write_bytes(b"accepted")
    (block / "trajectory.nc").write_bytes(b"trajectory")
    status = gates.read(block / "status.json")
    status.update(
        start_restart_sha256=gates.digest(previous),
        restart_sha256=gates.digest(block / "restart.ncrst"),
    )
    put(block / "status.json", status)
    process = gates.read(block / "process.json")
    process["start_restart_sha256"] = gates.digest(previous)
    put(block / "process.json", process)
    put(
        tmp_path / "extension_20ps_config.json",
        {"segment_input_sha256": {"extend_01": {}}},
    )
    return tmp_path


def test_actual_flat_native_receipt(native: Path) -> None:
    gates.accepted(native, 1)
    marker = native / "extend_01/inputs/restart.ncrst"
    marker.write_bytes(b"wrong predecessor")
    with pytest.raises(ValueError, match="provenance"):
        gates.accepted(native, 1)


@pytest.mark.parametrize(
    "field,value",
    [
        ("finite_temperature_passed", False),
        ("scc_verified", False),
        ("C_OG_A", float("nan")),
    ],
)
def test_reject_native_gates(native: Path, field: str, value: Any) -> None:
    path = native / "extend_01/scientific_acceptance.json"
    record = gates.read(path)
    record[field] = value
    put(path, record)
    with pytest.raises(ValueError):
        gates.accepted(native, 1)


def test_hashes_and_schema(tmp_path: Path) -> None:
    put(tmp_path / "list.json", [])
    with pytest.raises(ValueError):
        gates.read(tmp_path / "list.json")
    with pytest.raises(ValueError):
        gates.verify_hashes({})
    with pytest.raises(ValueError):
        gates.verify_hashes({str(tmp_path / "list.json"): "wrong"})


def test_initial_wait_and_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "extension_job/status.json"
    put(path, {"ended_utc": None})
    assert queue.initial_ready(tmp_path) is False
    put(path, {"ended_utc": "now", "exit_status": 1})
    with pytest.raises(ValueError):
        queue.initial_ready(tmp_path)
    put(path, {"ended_utc": "now", "exit_status": 0})
    put(
        tmp_path / "extension_20ps_status.json",
        {"status": "completed", "accepted_segments": 4},
    )
    monkeypatch.setattr(queue, "accepted", Mock())
    put(tmp_path / "analysis/report.json", {"frame_count_by_stage": {}})
    put(tmp_path / "refresh_recovery_job/status.json", {"ended_utc": None})
    assert queue.initial_ready(tmp_path) is False
    put(tmp_path / "analysis/report.json", {"frame_count_by_stage": {"extend_04": 100}})
    assert queue.initial_ready(tmp_path) is True


def test_exclusive_candidate_hold(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "analysis/source_manifest.json"
    put(path, {"sha256": {}})
    monkeypatch.setattr(queue, "report_candidates", lambda *_: [])
    assert queue.hold({"base": str(tmp_path)}, 4) is False
    monkeypatch.setattr(
        queue, "report_candidates", lambda *_: [{"near_attack_frame_count": 1}]
    )
    assert queue.hold({"base": str(tmp_path)}, 4) is True
    with pytest.raises(FileExistsError):
        queue.hold({"base": str(tmp_path)}, 4)


def test_report_actual_producer_schema(monkeypatch: pytest.MonkeyPatch) -> None:
    base = Path(
        "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
    ).resolve()
    # Real repaired report currently includes extend_01; native files retain their hashes.
    count = len(
        [
            s
            for s in gates.read(base / "analysis/report.json")["frame_count_by_stage"]
            if s.startswith("extend_")
        ]
    )
    assert isinstance(gates.report_candidates(base, count), list)
    with pytest.raises(ValueError, match="every accepted"):
        gates.report_candidates(base, count + 1)


def test_prepare_native_runtime(tmp_path: Path) -> None:
    original = Path(
        "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"
    )
    put(
        tmp_path / "extension_20ps_config.json",
        gates.read(original / "extension_20ps_config.json"),
    )
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for name in ("source.mdin", "system.prmtop"):
        (inputs / name).symlink_to((original / "inputs" / name).resolve())
    restart = tmp_path / "extend_04/restart.ncrst"
    restart.parent.mkdir()
    restart.write_bytes(b"accepted extend04 restart")
    snapshot = tmp_path / "followon_source_snapshot"
    snapshot.mkdir()
    (snapshot / "source.py").write_text("unchanged source")
    names = [f"extend_{i:02d}" for i in range(5, 10)]
    plan = {
        "base": str(tmp_path),
        "segment_names": names,
        "seeds": {name: 100 + i for i, name in enumerate(names)},
        "postprocess_argv": ["reviewed"],
        "postprocess_environment": {"PYTHONPATH": "explicit"},
        "recovery_source_sha256": {"reviewed": "source"},
    }
    queue.prepare_runtime(plan)
    cfg = gates.read(tmp_path / "followon_runtime_config.json")
    assert cfg["endpoint_restart_path"] == str(restart)
    assert (
        tmp_path / "extend_05/inputs/restart.ncrst"
    ).read_bytes() == restart.read_bytes()
    assert not (tmp_path / "extend_06/inputs/restart.ncrst").exists()
    native = (tmp_path / "extend_05/inputs/segment.mdin").read_text()
    assert "irest=1, ntx=5" in native and "nstlim=10000" in native
    with pytest.raises(FileExistsError):
        queue.prepare_runtime(plan)


@pytest.mark.parametrize(
    "mode",
    ["complete", "initial_candidate", "later_candidate", "native_fail", "timeout"],
)
def test_bounded_queue(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    from nyl12_covalent import deacylation_hydration_extension as native

    names = [f"extend_{i:02d}" for i in range(5, 10)]
    plan = {
        "base": str(tmp_path),
        "root": str(tmp_path),
        "wait_seconds": 0,
        "poll_seconds": 0,
        "segment_names": names,
        "pinned_sha256": {"key": "sha"},
    }
    path = tmp_path / "plan.json"
    put(path, plan)
    journal = Mock()
    monkeypatch.setattr(queue, "journal", journal)
    monkeypatch.setattr(queue, "verify_hashes", Mock())
    monkeypatch.setattr(queue, "initial_ready", lambda _: mode != "timeout")
    monkeypatch.setattr(
        queue,
        "hold",
        lambda _, n: (
            mode == "initial_candidate" or (mode == "later_candidate" and n == 5)
        ),
    )
    monkeypatch.setattr(queue.os, "cpu_count", lambda: 48)
    monkeypatch.setattr(queue.os, "getloadavg", lambda: (10.0, 10.0, 10.0))
    monkeypatch.setattr(
        queue.subprocess, "run", lambda *a, **k: Mock(stdout="resource snapshot")
    )
    monkeypatch.setattr(queue, "prepare_runtime", Mock())
    monkeypatch.setattr(queue, "accepted", Mock())
    run_segment = Mock(
        return_value={
            "status": "failed" if mode == "native_fail" else "accepted",
            "postprocess": {"status": "completed", "exit_code": 0},
        }
    )
    monkeypatch.setattr(native, "run_segment", run_segment)
    if mode in {"timeout", "native_fail"}:
        with pytest.raises((TimeoutError, ValueError)):
            queue.run(path)
    else:
        queue.run(path)
    assert (
        run_segment.call_count
        == {
            "complete": 5,
            "initial_candidate": 0,
            "later_candidate": 1,
            "native_fail": 1,
            "timeout": 0,
        }[mode]
    )
    assert journal.call_args.args[1] == "outcome"


def test_preparer_exclusive_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import shutil
    from scripts import prepare_hydration_followon as preparation

    original = preparation.BASE
    for name in ("extension_source_snapshot", "postprocess_recovery"):
        shutil.copytree(
            original / name,
            tmp_path / name,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
    for name in ("extension_20ps_config.json", "postprocess_recovery_config.json"):
        shutil.copyfile(original / name, tmp_path / name)
    monkeypatch.setattr(preparation, "BASE", tmp_path)
    preparation.prepare()
    plan = gates.read(tmp_path / "followon_25ps_plan.json")
    assert plan["total_limit_ps"] == 51 and len(plan["segment_names"]) == 5
    adapter = (
        tmp_path
        / "followon_postprocess_source/src/nyl12_covalent/hydration_receipt_adapter.py"
    )
    assert "range(1, 10)" in adapter.read_text()
    assert not (tmp_path / "extend_05").exists()
    with pytest.raises(FileExistsError):
        preparation.prepare()


def test_cli_entry(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from scripts import run_hydration_followon_queue as entry

    called = Mock()
    monkeypatch.setattr(entry, "run", called)
    path = tmp_path / "plan.json"
    monkeypatch.setattr(entry.sys, "argv", ["queue", "--plan", str(path)])
    entry.main()
    called.assert_called_once_with(path)
