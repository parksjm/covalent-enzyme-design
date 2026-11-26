"""Focused tests for reproducible hydration extension preparation and gates."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import deacylation_hydration as hydration
from nyl12_covalent import deacylation_hydration_extension as extension


def _pilot(tmp_path: Path) -> tuple[Path, dict[str, object]]:
    """Create a minimal hash-consistent accepted pilot fixture."""
    pilot = tmp_path / "hydration_reconnaissance"
    (pilot / "inputs").mkdir(parents=True)
    (pilot / "equil").mkdir()
    source = "&qmmm\nqmmask='@1-98', qmcharge=-2, qm_theory='DFTB3', dftb_telec=0.0, scfconv=1.d-12, qmshake=0,\n/\n&dlfind\nactive=':1-2',\n/\n"
    for name, payload in (
        ("endpoint.ncrst", b"restart"),
        ("system.prmtop", b"topology"),
        ("source.mdin", source.encode()),
    ):
        (pilot / "inputs" / name).write_bytes(payload)
    audit = tmp_path / "audit.json"
    audit.write_text("{}\n")
    executable = tmp_path / "sander"
    executable.write_text("executable")
    cfg = {
        "audit_path": str(audit),
        "audit_sha256": hydration.digest(audit),
        "input_sha256": {
            f"inputs/{name}": hydration.digest(pilot / "inputs" / name)
            for name in ("endpoint.ncrst", "system.prmtop", "source.mdin")
        },
        "sander": {"path": str(executable), "sha256": hydration.digest(executable)},
        "active_mask": ":1-2",
        "qmmask": "@1-98",
        "dftb_parameters_sha256": {},
        "journal_repo": str(tmp_path),
    }
    (pilot / "config.json").write_text(json.dumps(cfg))
    (pilot / "status.json").write_text(json.dumps({"status": "completed"}))
    acceptance = {
        "numeric": {
            "completed_steps": 10000,
            "final_temperature_K": 300.0,
            "electronic_history": {"verified_scc_history": True},
        },
        "boundary": {"frozen_max_displacement_A": 0.0, "box_unchanged": True},
        "chemistry": {"acyl_bond_retained": True, "C_OG_A": 1.31},
        "sampling_converged": False,
    }
    (pilot / "equil/process.json").write_text(
        json.dumps({"status": "completed", "exit_code": 0})
    )
    (pilot / "equil/scientific_acceptance.json").write_text(json.dumps(acceptance))
    endpoint_status = {
        "source_restart": "equil/restart.ncrst",
        "restart_sha256": hydration.digest(pilot / "inputs/endpoint.ncrst"),
        "pilot_status": {"status": "completed"},
        "final_acceptance": acceptance,
    }
    (pilot / "equil/restart.ncrst").write_bytes(
        (pilot / "inputs/endpoint.ncrst").read_bytes()
    )
    endpoint_status["restart_sha256"] = hydration.digest(pilot / "equil/restart.ncrst")
    (pilot / "hydration_endpoint_status_20261124.json").write_text(
        json.dumps(endpoint_status)
    )
    return pilot, cfg


def test_prepare_extension_preserves_pilot_and_builds_four_nvt_segments(
    tmp_path: Path,
) -> None:
    pilot, _ = _pilot(tmp_path)
    pilot_before = hydration.digest(pilot / "config.json")
    extension.prepare_extension(pilot, pilot, Path(__file__).resolve().parents[1])
    assert hydration.digest(pilot / "config.json") == pilot_before
    plan = json.loads((pilot / "extension_20ps_config.json").read_text())
    assert plan["segment_names"] == list(extension.SEGMENT_NAMES)
    assert plan["original_qmmm"] == {"atoms": 98, "charge": -2, "links": 8}
    assert plan["seed_provenance"].find("not claimed as bitwise") >= 0
    assert plan["postprocess_argv"] is None
    for name in extension.SEGMENT_NAMES:
        segment = pilot / name
        mdin = (segment / "inputs/segment.mdin").read_text()
        assert "nstlim=10000, dt=0.0005" in mdin
        assert "irest=1, ntx=5" in mdin
        assert "ntt=3, gamma_ln=2.0, temp0=300.0" in mdin
        assert "nmropt=0" in mdin and "ntr=0" in mdin
        assert "TEMP0" not in mdin and "DISANG" not in mdin
        assert f"ig={plan['seeds'][name]}" in mdin
        assert (segment / "inputs/system.prmtop").is_file()
        assert json.loads((segment / "status.json").read_text())["status"] == "prepared"
    assert hydration.digest(pilot / "extend_01/inputs/restart.ncrst") == plan["endpoint_restart_sha256"]
    assert not (pilot / "extend_02/inputs/restart.ncrst").exists()
    assert not (pilot / "heat").exists()
    with pytest.raises(FileExistsError):
        extension.prepare_extension(pilot, pilot, Path(__file__).resolve().parents[1])


def test_prepare_refuses_unaccepted_endpoint_before_staging(tmp_path: Path) -> None:
    pilot, _ = _pilot(tmp_path)
    endpoint_path = pilot / "hydration_endpoint_status_20261124.json"
    endpoint = json.loads(endpoint_path.read_text())
    endpoint["final_acceptance"]["chemistry"]["acyl_bond_retained"] = False
    endpoint_path.write_text(json.dumps(endpoint))
    with pytest.raises(ValueError, match="acceptance gates"):
        extension.prepare_extension(pilot, pilot, Path(__file__).resolve().parents[1])
    assert not (pilot / "extension_20ps_config.json").exists()


def test_input_hash_mismatch_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pilot, _ = _pilot(tmp_path)
    extension.prepare_extension(pilot, pilot, Path(__file__).resolve().parents[1])
    config = json.loads((pilot / "extension_20ps_config.json").read_text())
    monkeypatch.setattr(
        extension,
        "_load_identity",
        lambda path: (
            config,
            object(),
            np.ones(104922, dtype=bool),
            np.zeros((104922, 3)),
            np.ones(6),
        ),
    )
    target = pilot / "extend_01/inputs/segment.mdin"
    target.write_text(target.read_text() + "! changed\n")
    with pytest.raises(ValueError, match="input identity changed"):
        extension.run_segment(pilot, "extend_01")


def test_run_segment_accepts_native_and_scientific_gates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = tmp_path / "campaign"
    segment = campaign / "extend_01"
    (segment / "inputs").mkdir(parents=True)
    (segment / "status.json").write_text(
        json.dumps({"status": "prepared", "segment": 1})
    )
    (segment / "inputs/segment.mdin").write_text("mdin")
    (segment / "inputs/system.prmtop").write_text("top")
    endpoint = tmp_path / "endpoint.ncrst"
    endpoint.write_text("start")
    (segment / "inputs/restart.ncrst").write_text("start")
    cfg = {
        "segment_input_sha256": {
            "extend_01": {
                name: hydration.digest(segment / "inputs" / name)
                for name in ("segment.mdin", "system.prmtop", "restart.ncrst")
            }
        },
        "sander": {"path": "/bin/sander"},
        "seeds": {"extend_01": extension.SEEDS[0]},
        "postprocess_argv": None,
        "endpoint_restart_path": str(endpoint),
        "endpoint_restart_sha256": hydration.digest(endpoint),
        "journal_repo": str(tmp_path),
        "pilot_dir": str(tmp_path),
    }
    (campaign / "extension_20ps_config.json").write_text(json.dumps(cfg))
    active = np.zeros(104922, dtype=bool)
    active[:10] = True
    baseline = np.zeros((104922, 3))
    box = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    monkeypatch.setattr(
        extension, "_load_identity", lambda path: (cfg, object(), active, baseline, box)
    )
    monkeypatch.setattr(extension, "_verify_predecessor", lambda path, index: endpoint)
    monkeypatch.setattr(
        hydration,
        "audit_output",
        lambda text, steps: {
            "completed_steps": steps,
            "final_temperature_K": 300.0,
            "electronic_history": {"verified_scc_history": True},
        },
    )
    monkeypatch.setattr(
        hydration,
        "audit_restart",
        lambda *args: (
            np.zeros((104922, 3)),
            {"frozen_max_displacement_A": 0.0, "box_unchanged": True},
        ),
    )

    def fake_postprocess(
        campaign_path: Path, segment_path: Path, config: dict[str, object]
    ) -> None:
        (segment_path / "postprocess.json").write_text(
            json.dumps({"status": "completed"})
        )

    monkeypatch.setattr(extension, "_run_postprocess", fake_postprocess)

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        if command[0] == "/bin/sander":
            (segment / "run.mdout").write_text("native output")
            (segment / "restart.ncrst").write_text("new restart")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(extension.subprocess, "run", fake_run)
    result = extension.run_segment(campaign, "extend_01")
    assert result["status"] == "accepted"
    acceptance = json.loads((segment / "scientific_acceptance.json").read_text())
    assert acceptance["numeric"]["completed_steps"] == 10000
    assert acceptance["finite_temperature_passed"]
    assert acceptance["boundary"]["box_unchanged"]
    assert acceptance["acyl_retained"]
    assert (
        json.loads((segment / "postprocess.json").read_text())["status"] == "completed"
    )


def test_predecessor_must_be_accepted(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    (campaign / "extend_01").mkdir(parents=True)
    (campaign / "extend_01/status.json").write_text(json.dumps({"status": "prepared"}))
    (campaign / "extension_20ps_config.json").write_text(
        json.dumps({"endpoint_restart_path": str(tmp_path / "endpoint")})
    )
    with pytest.raises(ValueError, match="not scientifically accepted"):
        extension._verify_predecessor(campaign, 1)


def test_launch_refuses_unconfigured_postprocess_hook(tmp_path: Path) -> None:
    pilot, _ = _pilot(tmp_path)
    extension.prepare_extension(pilot, pilot, Path(__file__).resolve().parents[1])
    with pytest.raises(ValueError, match="source-pinned postprocess hook"):
        extension._load_identity(pilot)


def test_postprocess_records_hashes_and_failure_without_changing_md_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = tmp_path / "campaign"
    segment = campaign / "extend_01"
    (segment / "inputs").mkdir(parents=True)
    (campaign / "extension_source_snapshot/postprocess").mkdir(parents=True)
    source = tmp_path / "refresh.sh"
    source.write_text("#!/bin/sh\nexit 0\n")
    retained = campaign / "extension_source_snapshot/postprocess/refresh.sh"
    retained.write_bytes(source.read_bytes())
    source_hash = hydration.digest(source)
    config = {
        "journal_repo": str(Path(__file__).resolve().parents[1]),
        "postprocess_argv": ["/bin/bash", "{segment_dir}/refresh.sh", "{campaign}"],
        "postprocess_source_sha256": {str(source): source_hash},
    }
    (segment / "refresh.sh").write_bytes(source.read_bytes())
    calls: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        calls.append(command)
        if command[0] == "/bin/bash":
            source.write_text("changed after launch\n")
            return SimpleNamespace(returncode=0, stdout="export failed", stderr="")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(extension.subprocess, "run", fake_run)
    extension._run_postprocess(campaign, segment, config)
    receipt = json.loads((segment / "postprocess.json").read_text())
    assert receipt["status"] == "failed"
    assert "source identity changed" in receipt["error"]
    assert len(calls) == 3  # intent journal, configured hook, outcome journal
    assert hydration.digest(segment / "refresh.sh") == source_hash


def test_postprocess_success_records_source_hashes_and_exact_argv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = tmp_path / "campaign"
    segment = campaign / "extend_01"
    (segment / "inputs").mkdir(parents=True)
    retained_dir = campaign / "extension_source_snapshot/postprocess"
    retained_dir.mkdir(parents=True)
    source = tmp_path / "refresh.py"
    source.write_text("raise SystemExit(0)\n")
    retained = retained_dir / source.name
    retained.write_bytes(source.read_bytes())
    source_hash = hydration.digest(source)
    config = {
        "journal_repo": str(tmp_path),
        "postprocess_argv": ["python", str(source), "--base", "{campaign}"],
        "postprocess_source_sha256": {str(source): source_hash},
    }
    calls: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        calls.append(command)
        if command[0] == "python":
            return SimpleNamespace(returncode=0, stdout="refreshed", stderr="")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(extension.subprocess, "run", fake_run)
    extension._run_postprocess(campaign, segment, config)
    receipt = json.loads((segment / "postprocess.json").read_text())
    assert receipt["status"] == "completed"
    assert receipt["source_sha256_before"] == receipt["source_sha256_after"]
    assert receipt["argv"] == ["python", str(source), "--base", str(campaign.resolve())]
    assert len(calls) == 3


@pytest.mark.parametrize(
    ("mode", "expected"),
    [
        ("process", "process_failed"),
        ("scc", "scientific_failed"),
        ("acyl", "chemistry_stopped"),
    ],
)
def test_run_segment_stops_failed_process_or_scientific_gate(
    mode: str, expected: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = tmp_path / "campaign"
    segment = campaign / "extend_01"
    (segment / "inputs").mkdir(parents=True)
    for name, contents in (("segment.mdin", "mdin"), ("system.prmtop", "top")):
        (segment / "inputs" / name).write_text(contents)
    (segment / "status.json").write_text(json.dumps({"status": "prepared"}))
    endpoint = tmp_path / "endpoint.ncrst"
    endpoint.write_text("start")
    (segment / "inputs/restart.ncrst").write_text("start")
    cfg = {
        "segment_input_sha256": {
            "extend_01": {
                name: hydration.digest(segment / "inputs" / name)
                for name in ("segment.mdin", "system.prmtop")
            }
        },
        "sander": {"path": "/bin/sander"},
        "seeds": {"extend_01": extension.SEEDS[0]},
        "endpoint_restart_path": str(endpoint),
        "endpoint_restart_sha256": hydration.digest(endpoint),
        "journal_repo": str(tmp_path),
        "pilot_dir": str(tmp_path),
    }
    (campaign / "extension_20ps_config.json").write_text(json.dumps(cfg))
    active = np.ones(104922, dtype=bool)
    baseline = np.zeros((104922, 3))
    box = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    monkeypatch.setattr(
        extension, "_load_identity", lambda path: (cfg, object(), active, baseline, box)
    )
    monkeypatch.setattr(extension, "_verify_predecessor", lambda path, index: endpoint)

    def fake_output(text: str, steps: int) -> dict[str, object]:
        if mode == "scc":
            raise ValueError("SCC history failed")
        return {
            "completed_steps": steps,
            "final_temperature_K": 300.0,
            "electronic_history": {"verified_scc_history": True},
        }

    def fake_restart(*args: object) -> tuple[np.ndarray, dict[str, object]]:
        xyz = baseline.copy()
        if mode == "acyl":
            xyz[19915, 0] = 2.0
        return xyz, {"frozen_max_displacement_A": 0.0, "box_unchanged": True}

    monkeypatch.setattr(hydration, "audit_output", fake_output)
    monkeypatch.setattr(hydration, "audit_restart", fake_restart)
    monkeypatch.setattr(extension, "_run_postprocess", lambda *args: None)

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        if command[0] == "/bin/sander":
            if mode != "process":
                (segment / "run.mdout").write_text("native output")
                (segment / "restart.ncrst").write_text("restart")
            return SimpleNamespace(returncode=7 if mode == "process" else 0)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(extension.subprocess, "run", fake_run)
    result = extension.run_segment(campaign, "extend_01")
    assert result["status"] == expected
    assert json.loads((segment / "status.json").read_text())["status"] == expected


def test_extension_stops_on_first_failed_segment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[str] = []

    def run(campaign: Path, name: str) -> dict[str, object]:
        calls.append(name)
        return {"status": "accepted" if name == "extend_01" else "process_failed"}

    monkeypatch.setattr(extension, "run_segment", run)
    result = extension.run_extension(tmp_path)
    assert result["status"] == "stopped"
    assert calls == ["extend_01", "extend_02"]
    assert (
        json.loads((tmp_path / "extension_20ps_status.json").read_text())["segment"]
        == "extend_02"
    )


def test_extension_completes_only_after_all_four_accepted_segments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[str] = []

    def run(campaign: Path, name: str) -> dict[str, object]:
        calls.append(name)
        return {"status": "accepted", "segment": name}

    monkeypatch.setattr(extension, "run_segment", run)
    result = extension.run_extension(tmp_path)
    status = json.loads((tmp_path / "extension_20ps_status.json").read_text())
    assert result["status"] == status["status"] == "completed"
    assert status["accepted_segments"] == 4
    assert calls == list(extension.SEGMENT_NAMES)


def test_run_segment_rejects_unknown_or_already_started_segment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    with pytest.raises(ValueError, match="Unknown segment"):
        extension.run_segment(tmp_path, "segment_05")

    campaign = tmp_path / "campaign"
    segment = campaign / "extend_01"
    segment.mkdir(parents=True)
    (segment / "status.json").write_text(json.dumps({"status": "accepted"}))
    (campaign / "extension_20ps_config.json").write_text("{}")
    monkeypatch.setattr(
        extension,
        "_load_identity",
        lambda path: (
            {},
            object(),
            np.ones(2, dtype=bool),
            np.zeros((2, 3)),
            np.ones(6),
        ),
    )
    with pytest.raises(ValueError, match="no automatic retries"):
        extension.run_segment(campaign, "extend_01")
