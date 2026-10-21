"""Tests for reproducible relay-water repair job snapshots."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from scripts import launch_nyl12_relay_repair as launcher
from scripts import prepare_nyl12_inspection_download as export


def _repository(tmp_path: Path) -> tuple[Path, Path]:
    """Create the minimal repository layout required by the launcher."""
    root = tmp_path / "repo"
    scripts = root / "scripts"
    scripts.mkdir(parents=True)
    module_dir = root / "src/nyl12_covalent"
    module_dir.mkdir(parents=True)
    real_python = root / "real-python"
    real_python.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    real_python.chmod(0o755)
    python = root / ".venv/bin/python"
    python.parent.mkdir(parents=True)
    python.symlink_to(real_python)
    for name in ("repair_nyl12_relay_water.py", "launch_nyl12_relay_repair.py"):
        (scripts / name).write_text("# fixture\n", encoding="utf-8")
    for name in ("relay_repair.py", "relay_repair_geometry.py", "helper.py"):
        (module_dir / name).write_text("VALUE = 1\n", encoding="utf-8")
    (root / "pyproject.toml").write_text("[project]\n", encoding="utf-8")
    (root / "uv.lock").write_text("version = 1\n", encoding="utf-8")
    return root, python


def _capture_launch(
    monkeypatch: pytest.MonkeyPatch, captured: dict[str, Any]
) -> None:
    """Replace tmux integration with a deterministic captured call."""

    def fake_launch_job(
        session_name: str,
        output_dir: Path,
        cwd: Path,
        command: list[str],
        env_overrides: dict[str, str],
        input_paths: list[Path],
    ) -> Path:
        captured.update(
            session=session_name,
            output=output_dir,
            cwd=cwd,
            command=command,
            env=env_overrides,
            inputs=input_paths,
        )
        output_dir.mkdir(parents=True)
        return output_dir / "run.sh"

    monkeypatch.setattr(launcher, "launch_job", fake_launch_job)
    monkeypatch.setattr(
        launcher,
        "_runtime_versions",
        lambda python: {
            "python": "3.12.1",
            "executable": str(python),
            "numpy": "2.0.0",
            "parmed": "4.3.0",
        },
    )
    monkeypatch.setattr(
        launcher,
        "_git_snapshot",
        lambda repo: {
            "available": True,
            "head": "abc123",
            "status": " M tracked.py\n",
            "diff": b"tracked diff",
        },
    )


def test_launcher_snapshots_sources_inputs_and_uses_explicit_runtime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Job argv, environment, hashes, runtime, and source snapshots are explicit."""
    root, python = _repository(tmp_path)
    monkeypatch.setattr(launcher, "ROOT", root)
    topology = tmp_path / "same-name"
    coordinates = tmp_path / "other/same-name"
    coordinates.parent.mkdir()
    topology.write_bytes(b"topology")
    coordinates.write_bytes(b"coordinates")
    captured: dict[str, Any] = {}
    _capture_launch(monkeypatch, captured)

    attempt = tmp_path / "attempt"
    run_script = launcher.launch_repair_attempt(
        "pa6_l4", topology, coordinates, attempt, "relay_pa6_l4"
    )

    assert run_script == attempt / "job/run.sh"
    assert captured["output"] == attempt / "job"
    assert captured["cwd"] == attempt
    assert captured["session"] == "relay_pa6_l4"
    command = captured["command"]
    assert command[0] == str(python)
    assert command[1] == str(attempt / "source/scripts/repair_nyl12_relay_water.py")
    assert command[command.index("--target") + 1] == "pa6_l4"
    assert command[command.index("--step-degrees") + 1] == "15"
    assert command[command.index("--topology") + 1] == str(attempt / "inputs/ti1.prmtop")
    assert command[command.index("--coordinates") + 1] == str(attempt / "inputs/ti1.inpcrd")
    assert captured["env"] == {
        **launcher.THREAD_ENV,
        "PYTHONPATH": str(attempt / "source/src"),
    }
    manifest = json.loads((attempt / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["software_versions"]["numpy"] == "2.0.0"
    assert manifest["software_versions"]["parmed"] == "4.3.0"
    assert manifest["git"]["head"] == "abc123"
    assert (attempt / manifest["git"]["dirty_diff"]["path"]).read_bytes() == b"tracked diff"
    assert (attempt / manifest["git"]["status_snapshot"]["path"]).read_text() == (
        " M tracked.py\n"
    )
    snapshot_paths = {record["snapshot_relative_path"] for record in manifest["source_snapshot"]}
    assert "source/scripts/launch_nyl12_relay_repair.py" in snapshot_paths
    assert "source/pyproject.toml" in snapshot_paths
    assert "source/uv.lock" in snapshot_paths
    assert (attempt / "source/src/nyl12_covalent/helper.py").is_file()
    assert (attempt / "inputs/ti1.prmtop").read_bytes() == b"topology"
    assert (attempt / "inputs/ti1.inpcrd").read_bytes() == b"coordinates"
    for record in manifest["inputs"] + manifest["source_snapshot"]:
        original = Path(record["original_path"])
        snapshot = Path(record["snapshot_path"])
        assert hashlib.sha256(original.read_bytes()).hexdigest() == record["original_sha256"]
        assert hashlib.sha256(snapshot.read_bytes()).hexdigest() == record["snapshot_sha256"]


def test_preflight_failure_does_not_create_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Missing coordinates fail before a new attempt directory is made."""
    root, _ = _repository(tmp_path)
    monkeypatch.setattr(launcher, "ROOT", root)
    topology = tmp_path / "present.prmtop"
    topology.write_text("topology", encoding="utf-8")
    attempt = tmp_path / "attempt"
    with pytest.raises(FileNotFoundError):
        launcher.launch_repair_attempt(
            "maxa_maxa", topology, tmp_path / "missing.inpcrd", attempt, "relay_maxa"
        )
    assert not attempt.exists()


def test_runtime_version_and_git_helpers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Runtime and Git metadata come from their explicit subprocess results."""
    executable = tmp_path / "python"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    version_data = {
        "python": "3.12.1",
        "executable": str(executable),
        "numpy": "2.0.0",
        "parmed": "4.3.0",
    }
    monkeypatch.setattr(
        launcher.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout=json.dumps(version_data)),
    )
    assert launcher._runtime_versions(executable) == version_data

    outputs = [SimpleNamespace(stdout="head\n"), SimpleNamespace(stdout=" M x\n"), SimpleNamespace(stdout=b"diff")]
    monkeypatch.setattr(
        launcher.subprocess, "run", lambda *args, **kwargs: outputs.pop(0)
    )
    snapshot = launcher._git_snapshot(tmp_path)
    assert snapshot == {
        "available": True,
        "head": "head",
        "status": " M x\n",
        "diff": b"diff",
    }


def test_main_parses_required_arguments_and_reports_run_script(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The command-line entry point forwards all required launcher values."""
    received: dict[str, Any] = {}

    def fake_launch(
        target: str, topology: Path, coordinates: Path, output: Path, session: str
    ) -> Path:
        received.update(
            target=target,
            topology=topology,
            coordinates=coordinates,
            output=output,
            session=session,
        )
        return Path("attempt/job/run.sh")

    monkeypatch.setattr(launcher, "launch_repair_attempt", fake_launch)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "launcher",
            "--target",
            "maxa_maxa",
            "--topology",
            "top.prmtop",
            "--coordinates",
            "coords.inpcrd",
            "--output",
            "attempt",
            "--session",
            "relay_maxa",
        ],
    )
    assert launcher.main() == 0
    assert received == {
        "target": "maxa_maxa",
        "topology": Path("top.prmtop"),
        "coordinates": Path("coords.inpcrd"),
        "output": Path("attempt"),
        "session": "relay_maxa",
    }


def test_invalid_target_and_existing_attempt_fail_before_preflight(
    tmp_path: Path,
) -> None:
    """Invalid target and output collisions are rejected before input reads."""
    attempt = tmp_path / "attempt"
    with pytest.raises(ValueError, match="target must be one of"):
        launcher.launch_repair_attempt(
            "legacy", tmp_path / "missing.top", tmp_path / "missing.crd", attempt, "bad"
        )
    attempt.mkdir()
    with pytest.raises(ValueError, match="attempt directory already exists"):
        launcher.launch_repair_attempt(
            "pa6_l4", tmp_path / "missing.top", tmp_path / "missing.crd", attempt, "busy"
        )


def test_exporter_rejects_ambiguous_legacy_and_canonical_models(
    tmp_path: Path,
) -> None:
    """Both alias filenames cannot coexist in one source directory."""
    source = tmp_path / "models"
    source.mkdir()
    for name in (
        "pa6_l2_ti1_starting.pdb",
        "pa6_l4_ti1_starting.pdb",
        "maxa_maxa_ti1_starting.pdb",
    ):
        (source / name).write_text(name, encoding="utf-8")
    with pytest.raises(ValueError, match="Ambiguous source models for pa6_l4"):
        export.selected_files(
            tmp_path,
            tmp_path / "audit",
            tmp_path / "optimized",
            tmp_path / "trajectories",
            source,
        )


def test_exporter_can_package_a_frozen_action_journal(tmp_path: Path) -> None:
    """An explicit frozen journal takes precedence in the export inventory."""
    source = tmp_path / "models"
    source.mkdir()
    for name in ("pa6_l4_ti1_starting.pdb", "maxa_maxa_ti1_starting.pdb"):
        (source / name).write_text(name, encoding="utf-8")
    journal = tmp_path / "frozen-actions.jsonl"
    journal.write_text("{}\n", encoding="utf-8")
    inventory = export.selected_files(
        tmp_path,
        tmp_path / "audit",
        tmp_path / "optimized",
        tmp_path / "trajectories",
        source,
        journal,
    )
    assert inventory["provenance/agent_actions_at_export.jsonl"] == journal
