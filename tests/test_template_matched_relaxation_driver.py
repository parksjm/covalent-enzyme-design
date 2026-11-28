"""Checks for source-pinned ordering of the bounded native Amber controls."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
from types import ModuleType, SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot"
    / "wt_protein_template_controls_20261127/amber_relaxation_20261127"
    / "run_matched_relaxation.py"
)


@pytest.fixture
def driver() -> ModuleType:
    """Load frozen artifact without native execution."""
    spec = importlib.util.spec_from_file_location("matched_driver_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def plan(tmp_path: Path) -> dict:
    """Provide two new isolated case directories and a genuinely hashed input."""
    source = tmp_path / "source.pdb"
    source.write_text("immutable heavy input\n")
    return {
        "source_sha256": {str(source): hashlib.sha256(source.read_bytes()).hexdigest()},
        "config": {"gpu": "GPU-assigned"},
        "postprocessor": str(tmp_path / "postprocessor.py"),
        "cases": [
            {"name": name, "source": str(source), "directory": str(tmp_path / name)}
            for name in ("af3_template_seed20261126", "boltz_template_seed20261126")
        ],
    }


def test_source_change_rejected(driver: ModuleType, plan: dict) -> None:
    """A source edited after plan creation fails its retained SHA256 gate."""
    driver.verify_sources(plan)
    source = Path(next(iter(plan["source_sha256"])))
    source.write_text("changed source\n")
    with pytest.raises(ValueError, match="Source changed"):
        driver.verify_sources(plan)


@pytest.mark.parametrize(
    ("output", "blocked"),
    [
        ("", False),
        ("GPU-other, 123\n", False),
        ("GPU-other, 123\n GPU-assigned, 456\n", True),
    ],
)
def test_physical_gpu_clients(
    driver: ModuleType, monkeypatch: pytest.MonkeyPatch, output: str, blocked: bool
) -> None:
    """Other GPU clients are allowed; any assigned-GPU process blocks launch."""
    calls = []

    def query(command: list[str], **kwargs: object) -> SimpleNamespace:
        calls.append((command, kwargs))
        return SimpleNamespace(stdout=output)

    monkeypatch.setattr(driver.subprocess, "run", query)
    if blocked:
        with pytest.raises(RuntimeError, match="another compute client"):
            driver.require_idle_gpu("GPU-assigned")
    else:
        driver.require_idle_gpu("GPU-assigned")
    assert calls[0][0][0] == "nvidia-smi"
    assert calls[0][1]["check"] is True


@pytest.mark.parametrize("fail_second", [False, True])
def test_complete_native_pair_precedes_dry_exports(
    driver: ModuleType,
    plan: dict,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fail_second: bool,
) -> None:
    """Publish neither case if its matched native partner fails; route dry snapshots."""
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setattr(driver.sys, "argv", [str(SCRIPT), str(plan_path)])
    events = []
    exports = []

    def relax(source: Path, destination: Path, config: dict) -> dict:
        assert source == Path(plan["cases"][0]["source"])
        assert config == plan["config"]
        events.append(("relax", destination.parent.name))
        if fail_second and len(events) == 2:
            raise RuntimeError("native second case failed")
        return {"final_rms_gradient": 0.2, "force_converged": False}

    def execute(command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        if command[0] == "nvidia-smi":
            return subprocess.CompletedProcess(command, 0, stdout="GPU-other, 111\n")
        assert kwargs["check"] is True
        events.append(("export", command[command.index("--case") + 1]))
        exports.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(driver, "relax_model", relax)
    monkeypatch.setattr(driver.subprocess, "run", execute)
    if fail_second:
        with pytest.raises(RuntimeError, match="native second case failed"):
            driver.main()
        assert not exports
    else:
        driver.main()
        assert [event[0] for event in events] == ["relax", "relax", "export", "export"]
        for case, command in zip(plan["cases"], exports, strict=True):
            options = dict(zip(command[2::2], command[3::2], strict=True))
            native = Path(case["directory"]) / "relaxation"
            assert options == {
                "--case": case["name"],
                "--topology": str(native / "minimization/endpoint_complex.prmtop"),
                "--before-restart": str(native / "system/dry.rst7"),
                "--after-restart": str(native / "minimization/endpoint_complex.rst7"),
                "--conversion": str(native / "system/conversion.json"),
                "--output-dir": str(Path(case["directory"]) / "postprocess"),
            }
