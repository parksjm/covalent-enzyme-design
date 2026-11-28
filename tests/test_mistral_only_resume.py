"""Focused checks for the isolated Mistral recovery wrapper."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from scripts.run_mistral_only_resume import ISOLATED_PORT, run


def make_benchmark(
    root: Path, *, canary_ok: bool = True, full_ok: bool = True
) -> list[dict[str, Any]]:
    """Create a tiny fake frozen benchmark and capture model invocations."""
    source = root / "source_snapshot"
    source.mkdir(parents=True)
    (source / "protocol.py").write_text(
        "PORTS = {'nemotron': 18101, 'mistral': 18102}\n"
        "EXPECTED_MODELS = {'nemotron': {'port': 18101}, 'mistral': {'port': 18102}}\n"
        "def utc_now(): return 'test-time'\n",
        encoding="utf-8",
    )
    (root / "cases.json").write_text("{}", encoding="utf-8")
    settings = {"temperature": 0.7, "seed": 42}
    campaign = {
        "settings": settings,
        "system_prompt": "frozen prompt",
        "models": {
            "nemotron": {"model_name": "n", "port": 18101},
            "mistral": {"model_name": "m", "port": 18102},
        },
    }
    (root / "campaign.json").write_text(json.dumps(campaign), encoding="utf-8")
    calls: list[dict[str, Any]] = []
    runner_source = """
import json
from pathlib import Path

def load_inputs(root):
    campaign = json.loads((root / "campaign.json").read_bytes())
    return [{"case_id": "a"}, {"case_id": "b"}], campaign, "frozen-cases-hash"

def run_model(model, root, cases, campaign, cases_hash, resume):
    calls.append((model, [case["case_id"] for case in cases], campaign, cases_hash, resume))
    (root / "calls.json").write_text(json.dumps(calls))
    if len(calls) == 1 and not CANARY_OK:
        return False
    model_dir = root / "results" / "mistral"
    model_dir.mkdir(parents=True, exist_ok=True)
    entries = {}
    for case in cases:
        receipt = {"process_status": "complete", "finish_reason": "stop", "request_sha256": "hash-" + case["case_id"]}
        path = model_dir / (case["case_id"] + ".json")
        path.write_text(json.dumps(receipt))
        entries[case["case_id"]] = {"receipt": path.name}
    (model_dir / "progress.json").write_text(json.dumps({"cases": entries}))
    return True if len(cases) == 1 else FULL_OK
"""
    (source / "runner.py").write_text(
        "calls = CALLS\nCANARY_OK = CANARY_OK_VALUE\nFULL_OK = FULL_OK_VALUE\n"
        + runner_source,
        encoding="utf-8",
    )
    # Test constants are injected as literal Python values in this isolated fixture.
    text = (source / "runner.py").read_text(encoding="utf-8")
    text = (
        text.replace("CALLS", "[]")
        .replace("CANARY_OK_VALUE", repr(canary_ok))
        .replace("FULL_OK_VALUE", repr(full_ok))
    )
    (source / "runner.py").write_text(text, encoding="utf-8")
    return calls


def test_only_mistral_runs_with_transport_port_override(tmp_path: Path) -> None:
    """The wrapper retains scientific inputs and records isolated completion."""
    root = tmp_path / "benchmark"
    output = tmp_path / "recovery"
    make_benchmark(root)

    result = run(root, output)

    assert result["process"] == "complete"
    invocations = json.loads((root / "calls.json").read_text())
    assert [call[0:2] for call in invocations] == [
        ["mistral", ["a"]],
        ["mistral", ["a", "b"]],
    ]
    assert all(call[3:] == ["frozen-cases-hash", True] for call in invocations)
    assert all(
        call[2]["models"]["mistral"]["port"] == ISOLATED_PORT for call in invocations
    )
    assert all(
        call[2]["settings"] == {"temperature": 0.7, "seed": 42} for call in invocations
    )
    assert result["recovery_port"] == ISOLATED_PORT
    assert result["successful_stop_responses"] == 2
    assert result["request_sha256_by_case"] == {"a": "hash-a", "b": "hash-b"}
    assert (
        json.loads((output / "mistral_recovery_completion.json").read_text())["model"]
        == "mistral"
    )
    original = json.loads((root / "campaign.json").read_text())
    assert original["models"]["mistral"]["port"] == 18102
    assert original["settings"] == {"temperature": 0.7, "seed": 42}


def test_canary_failure_stops_before_full_resume(tmp_path: Path) -> None:
    """A failed first case prevents the full case list from being invoked."""
    root = tmp_path / "benchmark"
    make_benchmark(root, canary_ok=False)

    with pytest.raises(RuntimeError, match="canary"):
        run(root, tmp_path / "recovery")
    assert json.loads((root / "calls.json").read_text())[0][1] == ["a"]
    assert len(json.loads((root / "calls.json").read_text())) == 1


def test_rejects_any_other_endpoint_or_campaign_output(tmp_path: Path) -> None:
    """The recovery uses only its assigned port and an external receipt directory."""
    root = tmp_path / "benchmark"
    make_benchmark(root)

    with pytest.raises(ValueError, match="port"):
        run(root, tmp_path / "recovery", port=18102)
    with pytest.raises(ValueError, match="protected"):
        run(root, root / "results" / "recovery")


def test_new_recovery_directory_inside_campaign_is_allowed(tmp_path: Path) -> None:
    """A unique campaign-local recovery directory can retain its own receipt."""
    root = tmp_path / "benchmark"
    make_benchmark(root)
    output = root / "mistral_immediate_compat_retry01" / "benchmark_resume"

    result = run(root, output)

    assert result["process"] == "complete"
    assert (output / "mistral_recovery_completion.json").is_file()


def test_missing_progress_after_canary_fails_explicitly(tmp_path: Path) -> None:
    """The runner must publish progress after its canary returns successfully."""
    root = tmp_path / "benchmark"
    make_benchmark(root)
    runner_path = root / "source_snapshot" / "runner.py"
    source = runner_path.read_text(encoding="utf-8")
    source = source.replace(
        '    (model_dir / "progress.json").write_text(json.dumps({"cases": entries}))\n',
        "",
    )
    runner_path.write_text(source, encoding="utf-8")

    with pytest.raises(FileNotFoundError):
        run(root, tmp_path / "recovery")


def test_false_full_runner_result_records_failed_completion(tmp_path: Path) -> None:
    """A false runner result cannot be promoted to complete status."""
    root = tmp_path / "benchmark"
    make_benchmark(root, full_ok=False)

    result = run(root, tmp_path / "recovery")

    assert result["successful_stop_responses"] == 2
    assert result["process"] == "failed"
    saved = json.loads(
        (tmp_path / "recovery" / "mistral_recovery_completion.json").read_text()
    )
    assert saved["process"] == "failed"
