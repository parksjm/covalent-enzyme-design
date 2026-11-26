"""Verify no extension launches before reviewed, accepted baseline completion."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SPEC = importlib.util.spec_from_file_location(
    "extension_guard",
    Path(__file__).resolve().parents[1] / "scripts/guard_relay_umbrella_extension.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def fixture(tmp_path: Path) -> tuple[dict, dict, Path]:
    """Supply an explicit review and plan identity without starting Sander."""
    plan_path = tmp_path / "plan.json"
    plan_path.write_text("{}")
    plan = dict(
        original_plan=str(plan_path),
        original_plan_sha256=MODULE.digest(plan_path),
        sander="/test/sander",
        python="/test/python",
        workers=57,
        environment={},
    )
    config = dict(
        windows=[dict(label="w0")],
        topology="system.prmtop",
        temperature_K=300,
        physical_hamiltonian_sha256="a" * 64,
    )
    (tmp_path / "analysis_request.json").write_text(
        json.dumps(dict(diagnostic_columns=[]))
    )
    return plan, config, plan_path


def test_readiness_waits_for_future_outputs_and_rejects_failed_window(
    tmp_path: Path,
) -> None:
    plan, config, _ = fixture(tmp_path)
    result = MODULE.readiness(plan, config, tmp_path)
    assert result["ready"] is False and result["windows_ready"] == 0
    window = tmp_path / "w0"
    window.mkdir()
    (window / "window_status.json").write_text(json.dumps(dict(status="failed")))
    with pytest.raises(ValueError, match="window failed"):
        MODULE.readiness(plan, config, tmp_path)
    plan["original_plan_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="plan changed"):
        MODULE.readiness(plan, config, tmp_path)


def test_readiness_verifies_last_checkpoint_lineage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, config, _ = fixture(tmp_path)
    window = tmp_path / "w0"
    stage = window / "08_production_5ps"
    stage.mkdir(parents=True)
    previous = window / "07_production_5ps"
    previous.mkdir()
    (previous / "restart.ncrst").write_text("previous")
    (tmp_path / "system.prmtop").write_text("topology")
    for name in ["acceptance.json", "analysis_acceptance.json"]:
        (stage / name).write_text("{}")
    for name in ["restart.ncrst", "run.mdin", "umbrella.RST"]:
        (stage / name).write_text(name)
    (stage / "process.json").write_text(json.dumps(dict(exit_status=0)))
    inputs = [
        Path(plan["original_plan"]),
        tmp_path / "system.prmtop",
        previous / "restart.ncrst",
        stage / "run.mdin",
        stage / "umbrella.RST",
    ]
    (stage / "input_identity.json").write_text(
        json.dumps(dict(sha256={str(p): MODULE.digest(p) for p in inputs}))
    )
    (window / "segment_receipt_limit_06.json").write_text(
        json.dumps(
            dict(
                status="accepted",
                production_blocks=6,
                checkpoint_sha256=MODULE.digest(stage / "restart.ncrst"),
            )
        )
    )
    monkeypatch.setattr(MODULE, "load_block", lambda *a, **kw: dict(validated=True))
    assert MODULE.readiness(plan, config, tmp_path)["ready"] is True
    (stage / "restart.ncrst").write_text("changed")
    with pytest.raises(ValueError, match="restart changed"):
        MODULE.readiness(plan, config, tmp_path)


def test_launch_is_explicit_and_uses_correct_guard_cli(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, config, path = fixture(tmp_path)
    review = tmp_path / "review.json"
    review.write_text(json.dumps(dict(accepted=True, plan_sha256=MODULE.digest(path))))
    args = SimpleNamespace(plan=path, review=review, max_hours=72)
    monkeypatch.setattr(MODULE.argparse.ArgumentParser, "parse_args", lambda _: args)
    monkeypatch.setattr(MODULE, "verify_plan", lambda _: (plan, config, tmp_path))
    monkeypatch.setattr(
        MODULE, "readiness", lambda *a: dict(ready=True, windows_ready=57)
    )
    monkeypatch.setattr(MODULE, "native_active", lambda *a: [])
    calls = []
    monkeypatch.setattr(MODULE.subprocess, "run", lambda *a, **kw: calls.append(a[0]))
    jobs = []
    monkeypatch.setattr(MODULE, "launch_job", lambda *a, **kw: jobs.append(a))
    MODULE.main()
    assert len(jobs) == 1
    assert jobs[0][3][2] == "guard"
    assert "--analysis-only" in jobs[0][3]
    assert "--poll-seconds" not in jobs[0][3]
    assert (
        json.loads((tmp_path / "gate_status.json").read_text())["native_launched"]
        is True
    )
    (tmp_path / "launch_claim").write_text("claimed")
    with pytest.raises(FileExistsError):
        MODULE.main()
    review.write_text(json.dumps(dict(accepted=False)))
    with pytest.raises(ValueError, match="reviewed"):
        MODULE.main()
    args.max_hours = 0
    with pytest.raises(ValueError, match="Positive"):
        MODULE.main()
