"""Continuation requires a clean closed original queue and finalized claims."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
from unittest.mock import Mock

import pytest
from scripts import guard_relay_umbrella as guardian
from scripts.guard_relay_umbrella import inspect_held


def test_guard_only_continues_clean_held_queue(tmp_path: Path) -> None:
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps(dict(workers=1)))
    folder = tmp_path / "workers/worker_00"
    (folder / "job").mkdir(parents=True)
    status = folder / "job/status.json"
    status.write_text(json.dumps(dict(ended_utc=None, exit_status=None)))
    assert inspect_held(plan) is False
    status.write_text(json.dumps(dict(ended_utc="closed", exit_status=0)))
    (folder / "status.json").write_text(json.dumps(dict(status="held_for_analysis")))
    (tmp_path / "claims").mkdir()
    claim = tmp_path / "claims/a.json"
    claim.write_text(json.dumps(dict(status="completed")))
    assert inspect_held(plan) is True
    claim.write_text(json.dumps(dict(status="running")))
    with pytest.raises(RuntimeError, match="Unresolved native queue claim"):
        inspect_held(plan)
    claim.write_text(json.dumps(dict(status="completed")))
    status.write_text(json.dumps(dict(ended_utc="closed", exit_status=1)))
    with pytest.raises(RuntimeError, match="Native worker failure"):
        inspect_held(plan)


def test_analysis_only_repairs_receipts_without_touching_native_queue(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import parmed
    from nyl12_covalent import jobs, relay_umbrella_analysis_io, relay_umbrella_runtime

    stage = tmp_path / "window/04_production_5ps"
    stage.mkdir(parents=True)
    native = stage / "acceptance.json"
    native.write_text('{"accepted": true}\n')
    original_native = native.read_bytes()
    (stage / "analysis_outcome.json").write_text('{"status": "failed"}\n')
    (stage.parent / "segment_receipt_limit_01.json").write_text("{}\n")
    (tmp_path / "prepared.json").write_text(
        json.dumps(
            dict(
                topology="system.prmtop",
                temperature_K=300,
                physical_hamiltonian_sha256="unchanged",
                windows=[
                    dict(
                        label="window",
                        stages=[
                            dict(
                                kind="production",
                                label="04_production_5ps",
                                input="window/04_production_5ps/run.mdin",
                            )
                        ],
                    )
                ],
            )
        )
    )
    (tmp_path / "analysis_request.json").write_text('{"diagnostic_columns": []}\n')
    plan = tmp_path / "plan.json"
    plan.write_text('{"workers": 57}\n')
    calls = iter((0.0, 0.0, 2.0))
    monkeypatch.setattr(guardian.time, "monotonic", lambda: next(calls))
    monkeypatch.setattr(guardian.time, "sleep", lambda _: None)
    monkeypatch.setattr(guardian, "journal", lambda *_: None)
    native_queue = Mock(side_effect=AssertionError("Native queue must be untouched"))
    monkeypatch.setattr(guardian, "inspect_held", native_queue)
    monkeypatch.setattr(jobs, "launch_job", native_queue)
    topology = object()
    monkeypatch.setattr(parmed, "load_file", lambda _: topology)
    repair = Mock(return_value=True)
    monkeypatch.setattr(relay_umbrella_runtime, "ensure_diagnostics", repair)
    accepted = Mock()
    monkeypatch.setattr(relay_umbrella_analysis_io, "load_block", accepted)
    output = tmp_path / "repair_state"

    guardian.guard(
        argparse.Namespace(
            campaign=tmp_path,
            plan=plan,
            output=output,
            duration_hours=1 / 3600,
            analysis_only=True,
        )
    )

    repair.assert_called_once()
    assert repair.call_args.args[2:] == (stage, topology)
    accepted.assert_called_once()
    native_queue.assert_not_called()
    assert native.read_bytes() == original_native
    assert (
        json.loads((output / "latest.json").read_text())["repairs"][0][
            "analysis_completed"
        ]
        is True
    )


def test_analysis_only_is_rejected_for_native_worker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "guard_relay_umbrella.py",
            "worker",
            "--plan",
            str(tmp_path / "plan.json"),
            "--output",
            str(tmp_path / "output"),
            "--index",
            "0",
            "--analysis-only",
        ],
    )
    with pytest.raises(SystemExit, match="2"):
        guardian.main()
