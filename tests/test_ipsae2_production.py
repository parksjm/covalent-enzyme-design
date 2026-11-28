"""Production preparation and terminal watcher contract tests."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from nyl12_covalent import ipsae2_prepare as prepare
from nyl12_covalent import ipsae2_scoring as scoring
from nyl12_covalent import ipsae2_sidecar as sidecar


def write(path: Path, data: Any) -> Path:
    """Write a retained synthetic JSON input."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return path


@pytest.fixture
def prepared_inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """Build six complete native directories and mock only external tools."""
    executable = tmp_path / "ipsae"
    executable.write_text("pinned")
    runtime = {
        "command": [str(executable)],
        "revision": scoring.REVISION,
        "version": "0.1.11",
        "source_sha256": {str(executable): scoring.file_sha256(executable)},
    }
    runtime_path = write(tmp_path / "runtime.json", runtime)
    qualification = write(
        tmp_path / "qualification.json",
        {
            "status": "accepted",
            "runtime_sha256": scoring.runtime_digest(runtime),
            "qualified_roles": ["WT", "design"],
        },
    )
    campaign = tmp_path / "campaign"
    for index in range(6):
        case = f"WT_seed_{index}"
        base = campaign / "prediction/cases" / case
        result = base / "runs" / f"boltz_results_{case}"
        prediction = result / "predictions" / case
        prediction.mkdir(parents=True)
        cif = prediction / f"{case}_model_0.cif"
        cif.write_text("native")
        for kind, suffix in (("pae", "npz"), ("plddt", "npz"), ("confidence", "json")):
            (prediction / f"{kind}_{cif.stem}.{suffix}").write_text(kind)
        structure = result / f"processed/structures/{case}.npz"
        structure.parent.mkdir(parents=True)
        structure.write_text("structure")
        write(
            base / "native_prediction.json",
            {
                "status": "native_prediction_completed",
                "case": case,
                "candidate": "WT",
                "seed": index,
                "prediction_cif": str(cif),
                "prediction_cif_sha256": scoring.file_sha256(cif),
            },
        )
    for file in (
        "src/nyl12_covalent/ipsae2_scoring.py",
        "src/nyl12_covalent/ipsae2_sidecar.py",
        "scripts/score_boltz_ipsae2.py",
    ):
        path = tmp_path / file
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("source")

    def mapping_run(command: list[str], **kwargs: Any) -> SimpleNamespace:
        artifacts = {}
        for key in ("cif", "pae", "plddt", "confidence", "structure"):
            path = command[command.index(f"--{key}") + 1]
            artifacts[key] = {"path": path, "sha256": scoring.file_sha256(path)}
        write(
            Path(command[command.index("--output") + 1]),
            {"accepted": True, "artifacts": artifacts},
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(prepare.subprocess, "run", mapping_run)
    monkeypatch.setattr(
        prepare,
        "physical_pairing",
        lambda *a: {
            "protomer_groups": [[a, b] for a, b in zip("ABCD", "EFGH", strict=True)]
        },
    )
    return {
        "campaign": campaign,
        "output": tmp_path / "scoring",
        "runtime_path": runtime_path,
        "qualification_path": qualification,
        "pairing_adapter": tmp_path / "adapter.py",
        "boltz_python": tmp_path / "boltz_python",
        "boltz_source_root": tmp_path / "boltz",
        "root": tmp_path,
    }


def test_prepare_all_six(prepared_inputs: dict[str, Path]) -> None:
    config_path = prepare.prepare_completed_campaign(**prepared_inputs)
    config = sidecar.read_json(config_path)
    assert len(config["cases"]) == 6
    assert all(Path(case["mapping_receipt"]).exists() for case in config["cases"])
    assert all(Path(case["model_manifest"]).exists() for case in config["cases"])
    assert len(config["input_sha256"]) == 29
    assert all(
        scoring.file_sha256(path) == value
        for path, value in config["input_sha256"].items()
    )
    assert not (prepared_inputs["output"] / "WT_seed_0").exists()


@pytest.mark.parametrize("failure", ["count", "pending", "identity", "hash", "mapping"])
def test_prepare_rejects_invalid(
    prepared_inputs: dict[str, Path], failure: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    native = next(
        (prepared_inputs["campaign"] / "prediction/cases").glob(
            "*/native_prediction.json"
        )
    )
    record = sidecar.read_json(native)
    if failure == "count":
        native.unlink()
    elif failure == "pending":
        record["status"] = "running"
        write(native, record)
    elif failure == "identity":
        record["case"] = "different"
        write(native, record)
    elif failure == "hash":
        Path(record["prediction_cif"]).write_text("tampered")
    else:
        monkeypatch.setattr(
            prepare.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=1)
        )
    with pytest.raises((ValueError, RuntimeError)):
        prepare.prepare_completed_campaign(**prepared_inputs)


def test_native_failure_terminal_and_stale_pin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sidecar, "verify_runtime", lambda *a: None)
    monkeypatch.setattr(sidecar, "require_qualification", lambda *a: None)
    path = write(
        tmp_path / "native.json",
        {"status": "native_prediction_failed", "exit_status": 1},
    )
    config = {
        "runtime": str(write(tmp_path / "runtime.json", {})),
        "qualification": str(write(tmp_path / "q.json", {})),
        "scoring_root": str(tmp_path / "scores"),
        "cases": [{"case_id": "WT_1", "native_receipt": str(path)}],
    }
    assert sidecar.poll_cases(config)["cases"]["WT_1"]["status"] == "failed"
    config["input_sha256"] = {str(path): "0" * 64}
    with pytest.raises(ValueError, match="Sidecar input hash"):
        sidecar.poll_cases(config)


def test_watch_exits_when_every_case_terminal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    script = Path(__file__).resolve().parents[1] / "scripts/score_boltz_ipsae2.py"
    specification = importlib.util.spec_from_file_location("score_boltz_ipsae2", script)
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    monkeypatch.setattr(
        sys, "argv", [str(script), "--config", str(tmp_path / "config"), "--watch"]
    )
    monkeypatch.setattr(module, "read_json", lambda *a: {})
    monkeypatch.setattr(
        module,
        "poll_cases",
        lambda *a, **k: {
            "cases": {"WT": {"status": "completed"}, "design": {"status": "failed"}}
        },
    )
    monkeypatch.setattr(
        module.time, "sleep", lambda *a: pytest.fail("terminal watcher slept")
    )
    module.main()
