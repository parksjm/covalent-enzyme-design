"""Exercise conditional w033 gate and watcher without launching science."""

import json

import pytest

from nyl12_covalent.relay_path_queue import digest, verify_minimum
from scripts import watch_window033_relay as watcher


@pytest.fixture
def evidence(tmp_path):
    frozen = tmp_path / "01_frozen"
    files = (
        "inputs/optimized.ncrst",
        "inputs/system.prmtop",
        "inputs/dimer.mdin",
        "library_build/libdlfind.so",
    )
    for name in files:
        path = frozen / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    (frozen / "manifest.json").write_text(
        json.dumps({n: digest(frozen / n) for n in files})
    )
    stage = tmp_path / "06_full_hessian"
    (stage / "job").mkdir(parents=True)
    (stage / "job/status.json").write_text(
        json.dumps({"ended_utc": "done", "exit_status": 0})
    )
    (stage / "specification.json").write_text(
        json.dumps(
            {
                "dimension": 12495,
                "restart_sha256": digest(frozen / files[0]),
                "library_sha256": digest(frozen / files[3]),
            }
        )
    )
    inertia = {
        "passed": True,
        "negative": 0,
        "unresolved": 0,
        "expected_negative": 0,
        "dimension": 12495,
        "raw_hessian_sha256": "raw",
        "symmetric_hessian_sha256": "symmetric",
    }
    (stage / "scientific_acceptance.json").write_text(json.dumps(inertia))
    (tmp_path / "verification_summary.json").write_text(
        json.dumps(
            {
                "status": "numerical_local_minimum",
                "numerically_verified": True,
                "literal_inertia": inertia,
            }
        )
    )
    return tmp_path


def update(path, **changes):
    data = json.loads(path.read_text())
    data.update(changes)
    path.write_text(json.dumps(data))


def test_matching_receipt(evidence):
    result = verify_minimum(evidence)
    assert result["restart_sha256"] == digest(
        evidence / "01_frozen/inputs/optimized.ncrst"
    )
    assert result["hessian_acceptance_sha256"] == digest(
        evidence / "06_full_hessian/scientific_acceptance.json"
    )
    assert result["literal_inertia"]["dimension"] == 12495


def test_pending_until_terminal_process(evidence):
    status = evidence / "06_full_hessian/job/status.json"
    update(status, ended_utc=None, exit_status=None)
    assert verify_minimum(evidence) is None
    update(
        evidence / "verification_summary.json",
        status="full_hessian_running",
        numerically_verified=False,
    )
    assert verify_minimum(evidence) is None


@pytest.mark.parametrize(
    "kind",
    [
        "not_verified",
        "incomplete",
        "failed_process",
        "negative",
        "unresolved",
        "expected",
        "dimension",
        "mismatch",
        "hash",
        "spec",
    ],
)
def test_rejects_invalid_evidence(evidence, kind):
    summary = evidence / "verification_summary.json"
    report = json.loads(summary.read_text())
    if kind in ("not_verified", "incomplete"):
        report["status"] = kind
    elif kind == "failed_process":
        update(evidence / "06_full_hessian/job/status.json", exit_status=1)
    elif kind in ("negative", "unresolved", "expected", "dimension"):
        key = {"expected": "expected_negative"}.get(kind, kind)
        report["literal_inertia"][key] = 1
    elif kind == "mismatch":
        update(evidence / "06_full_hessian/scientific_acceptance.json", dimension=3)
    elif kind == "hash":
        (evidence / "01_frozen/inputs/optimized.ncrst").write_text("changed")
    else:
        update(
            evidence / "06_full_hessian/specification.json", restart_sha256="changed"
        )
    summary.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        verify_minimum(evidence)


def test_watcher_dispatch_and_reuse_rejected(evidence, tmp_path, monkeypatch):
    campaign = tmp_path / "campaign"
    driver = campaign / "source_snapshot/scripts/driver.py"
    driver.parent.mkdir(parents=True)
    driver.write_text("pass")
    (campaign / "config.json").write_text(
        json.dumps({"verification_root": str(evidence)})
    )
    calls = []
    monkeypatch.setattr(watcher, "_journal", lambda *args: None)
    monkeypatch.setattr(
        watcher.subprocess, "call", lambda argv: calls.append(argv) or 0
    )
    assert watcher.watch(evidence, campaign, driver, 1, 1) == 0
    assert calls[0][1:] == [str(driver), "--campaign", str(campaign)]
    assert (
        json.loads((campaign / "minimum_gate_watcher/wait.json").read_text())["status"]
        == "completed"
    )
    assert [
        json.loads(row)["event"]
        for row in (campaign / "minimum_gate_watcher/events.jsonl")
        .read_text()
        .splitlines()
    ] == ["watch_started", "dispatch_intent", "dispatch_outcome"]
    with pytest.raises(FileExistsError):
        watcher.watch(evidence, campaign, driver, 1, 1)


def test_watcher_timeout_and_mutation(evidence, tmp_path, monkeypatch):
    campaign = tmp_path / "pending"
    driver = campaign / "source_snapshot/driver.py"
    driver.parent.mkdir(parents=True)
    driver.write_text("pass")
    config = campaign / "config.json"
    config.write_text(json.dumps({"verification_root": str(evidence)}))
    monkeypatch.setattr(watcher, "verify_minimum", lambda root: None)
    monkeypatch.setattr(
        watcher.time, "sleep", lambda value: driver.write_text("changed")
    )
    with pytest.raises(ValueError, match="changed"):
        watcher.watch(evidence, campaign, driver, 0.01, 1)
    assert (
        json.loads((campaign / "minimum_gate_watcher/wait.json").read_text())["status"]
        == "failed"
    )


def test_watcher_invalid_inputs(evidence, tmp_path):
    with pytest.raises(ValueError, match="positive"):
        watcher.watch(evidence, tmp_path, tmp_path, 0, 1)
    campaign = tmp_path / "bad"
    (campaign / "source_snapshot").mkdir(parents=True)
    driver = tmp_path / "outside.py"
    driver.write_text("pass")
    with pytest.raises(ValueError, match="source snapshot"):
        watcher.watch(evidence, campaign, driver, 1, 1)


def test_gate_process_failure_without_summary(evidence):
    (evidence / "verification_summary.json").unlink()
    update(evidence / "06_full_hessian/job/status.json", exit_status=3)
    with pytest.raises(ValueError, match="native job failed"):
        verify_minimum(evidence)


def test_journal_failure_blocks_dispatch(evidence, tmp_path, monkeypatch):
    campaign = tmp_path / "journal"
    driver = campaign / "source_snapshot/driver.py"
    driver.parent.mkdir(parents=True)
    driver.write_text("pass")
    (campaign / "config.json").write_text(
        json.dumps({"verification_root": str(evidence)})
    )
    calls = []
    monkeypatch.setattr(watcher.subprocess, "call", lambda argv: calls.append(argv))

    def reject(*args):
        raise RuntimeError("journal unavailable")

    monkeypatch.setattr(watcher, "_journal", reject)
    with pytest.raises(RuntimeError, match="journal unavailable"):
        watcher.watch(evidence, campaign, driver, 1, 1)
    assert calls == []
