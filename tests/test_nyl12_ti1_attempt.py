from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from scripts import build_nyl12_ti1_qmmm
from nyl12_covalent import qmmm, ti1_attempt


@pytest.fixture
def attempt_inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """Create small immutable inputs and a synthetic Amber installation."""
    topology = tmp_path / "source.prmtop"
    coordinates = tmp_path / "source.rst7"
    reference = tmp_path / "reference.rst7"
    for path, content in (
        (topology, b"topology-v1\n"),
        (coordinates, b"coordinates-v1\n"),
        (reference, b"reference-v1\n"),
    ):
        path.write_bytes(content)

    amber_bin = tmp_path / "amber" / "bin"
    amber_bin.mkdir(parents=True)
    sander = amber_bin / "sander"
    sander.write_text("#!/bin/sh\nexit 0\n")
    sander.chmod(sander.stat().st_mode | 0o111)
    slko = tmp_path / "amber" / "dat" / "slko"
    slko.mkdir(parents=True)
    (slko / "parameters.skf").write_text("synthetic DFTB parameter\n")

    monkeypatch.setattr(ti1_attempt, "_snapshot", lambda _output: {"test.py": "abc"})
    monkeypatch.setattr(
        ti1_attempt,
        "audit_qm_region",
        lambda _topology: {"accepted": True, "qmmask": qmmm.QM_MASK},
    )
    monkeypatch.setattr(
        qmmm,
        "audit_qm_region",
        lambda _topology: {"accepted": True, "qmmask": qmmm.QM_MASK},
    )
    monkeypatch.setattr(
        ti1_attempt,
        "check_geometry",
        lambda *_args: {"accepted": True, "screen": "pass"},
    )
    return {
        "topology": topology,
        "coordinates": coordinates,
        "reference": reference,
        "amber_bin": amber_bin,
        "sander": sander,
        "slko": slko,
        "root": tmp_path,
    }


def _run(
    inputs: dict[str, Path],
    output: Path,
    *,
    coordinates: Path | None = None,
    checkpoint_attempt: Path | None = None,
    target: str = "pa6_l2",
) -> Path:
    return ti1_attempt.run_attempt(
        topology=inputs["topology"],
        coordinates=coordinates or inputs["coordinates"],
        reference=inputs["reference"],
        output=output,
        amber_bin=inputs["amber_bin"],
        target=target,
        checkpoint_attempt=checkpoint_attempt,
    )


def _mock_engine(
    monkeypatch: pytest.MonkeyPatch,
    *,
    returncode: int = 0,
    mdout: str = "engine output\n",
) -> list[list[str]]:
    commands: list[list[str]] = []

    def run(
        command: list[str],
        *,
        cwd: Path,
        stdout: Any,
        **_kwargs: Any,
    ) -> SimpleNamespace:
        commands.append(command)
        Path(cwd, "ti1_dftb3.rst7").write_bytes(b"engine-restart\n")
        Path(cwd, "ti1_dftb3.mdout").write_text(mdout)
        stdout.write("synthetic sander log\n")
        return SimpleNamespace(returncode=returncode)

    monkeypatch.setattr(ti1_attempt.subprocess, "run", run)
    return commands


def test_failed_preflight_preserves_inputs_and_does_not_start_engine(
    attempt_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    original = {
        name: attempt_inputs[name].read_bytes()
        for name in ("topology", "coordinates", "reference")
    }
    monkeypatch.setattr(
        ti1_attempt,
        "check_geometry",
        lambda *_args: {"accepted": False, "reason": "relay screen failed"},
    )
    engine = _mock_engine(monkeypatch)
    output = attempt_inputs["root"] / "failed-preflight"

    with pytest.raises(ValueError, match="Starting geometry failed"):
        _run(attempt_inputs, output)

    assert engine == []
    assert {name: attempt_inputs[name].read_bytes() for name in original} == original
    assert (output / "ti1.prmtop").read_bytes() == original["topology"]
    assert json.loads((output / "acceptance.json").read_text())["stage"] == "preflight"
    process = json.loads((output / "process.json").read_text())
    assert process["ended_utc"] is not None
    assert process["exit_status"] is None


def test_process_without_convergence_is_rejected_and_recorded(
    attempt_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    engine = _mock_engine(monkeypatch, mdout="no final results\n")
    output = attempt_inputs["root"] / "no-convergence"

    with pytest.raises(ValueError, match="convergence/chemistry acceptance failed"):
        _run(attempt_inputs, output)

    assert len(engine) == 1
    acceptance = json.loads((output / "acceptance.json").read_text())
    assert acceptance["convergence"]["accepted"] is False
    assert acceptance["accepted"] is False
    assert json.loads((output / "process.json").read_text())["exit_status"] == 0


def test_fresh_output_directory_is_never_overwritten(
    attempt_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    engine = _mock_engine(monkeypatch)
    output = attempt_inputs["root"] / "existing"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("keep me")

    with pytest.raises(FileExistsError):
        _run(attempt_inputs, output)

    assert sentinel.read_text() == "keep me"
    assert engine == []


def _successful_attempt(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch, output: Path
) -> Path:
    _mock_engine(monkeypatch)
    monkeypatch.setattr(
        ti1_attempt,
        "check_convergence",
        lambda _mdout: {"accepted": True, "final_results_present": True},
    )
    return _run(inputs, output)


@pytest.mark.parametrize(
    "drift", ["topology", "reference", "specification", "checkpoint"]
)
def test_continuation_rejects_scientific_or_checkpoint_drift(
    attempt_inputs: dict[str, Path],
    monkeypatch: pytest.MonkeyPatch,
    drift: str,
) -> None:
    source = attempt_inputs["root"] / "source-attempt"
    _successful_attempt(attempt_inputs, monkeypatch, source)
    checkpoint = source / "ti1_dftb3.rst7"
    target = "pa6_l2"
    if drift == "topology":
        attempt_inputs["topology"].write_text("topology-v2\n")
    elif drift == "reference":
        attempt_inputs["reference"].write_text("reference-v2\n")
    elif drift == "specification":
        manifest_path = source / "attempt.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["specification"]["qmcharge"] = 1
        manifest_path.write_text(json.dumps(manifest))
    else:
        checkpoint.write_text("modified checkpoint\n")
    output = attempt_inputs["root"] / f"continued-{drift}"

    with pytest.raises(ValueError, match="mismatch|corrupt|hash"):
        _run(
            attempt_inputs,
            output,
            coordinates=checkpoint,
            checkpoint_attempt=source,
            target=target,
        )

    assert not output.exists()


def test_accepted_attempt_returns_output_and_records_checkpoint_hashes(
    attempt_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    output = attempt_inputs["root"] / "accepted"
    mdout = _successful_attempt(attempt_inputs, monkeypatch, output)

    assert mdout == output / "ti1_dftb3.mdout"
    assert mdout.is_file()
    manifest = json.loads((output / "attempt.json").read_text())
    process = json.loads((output / "process.json").read_text())
    acceptance = json.loads((output / "acceptance.json").read_text())
    assert manifest["inputs"]["topology"]["sha256"]
    assert process["restart_sha256"] == ti1_attempt._sha(output / "ti1_dftb3.rst7")
    assert acceptance["accepted"] is True


def test_successful_continuation_uses_checkpoint_and_original_reference(
    attempt_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    source = attempt_inputs["root"] / "source-attempt"
    _successful_attempt(attempt_inputs, monkeypatch, source)
    source_restart = source / "ti1_dftb3.rst7"
    source_reference = source / "reference.rst7"
    source_reference_hash = ti1_attempt._sha(source_reference)
    commands = _mock_engine(monkeypatch)
    monkeypatch.setattr(
        ti1_attempt,
        "check_convergence",
        lambda _mdout: {"accepted": True, "final_results_present": True},
    )
    continued = attempt_inputs["root"] / "continued"

    result = _run(
        attempt_inputs,
        continued,
        coordinates=source_restart,
        checkpoint_attempt=source,
    )

    assert result == continued / "ti1_dftb3.mdout"
    command = commands[0]
    start_index = command.index("-c") + 1
    reference_index = command.index("-ref") + 1
    assert Path(command[start_index]) == continued / "start.rst7"
    assert ti1_attempt._sha(Path(command[start_index])) == ti1_attempt._sha(
        source_restart
    )
    assert Path(command[reference_index]) == continued / "reference.rst7"
    assert ti1_attempt._sha(Path(command[reference_index])) == source_reference_hash
    provenance = json.loads((continued / "attempt.json").read_text())
    assert provenance["checkpoint_source"]["checkpoint_sha256"] == ti1_attempt._sha(
        source_restart
    )


def test_continuation_rejects_source_without_terminal_process_record(
    attempt_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    source = attempt_inputs["root"] / "unfinished-source"
    _successful_attempt(attempt_inputs, monkeypatch, source)
    process_path = source / "process.json"
    process = json.loads(process_path.read_text())
    process["ended_utc"] = None
    process_path.write_text(json.dumps(process))
    output = attempt_inputs["root"] / "must-not-start"

    with pytest.raises(ValueError, match="no terminal process record"):
        _run(
            attempt_inputs,
            output,
            coordinates=source / "ti1_dftb3.rst7",
            checkpoint_attempt=source,
        )

    assert not output.exists()


def test_engine_failure_writes_terminal_records(
    attempt_inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    _mock_engine(monkeypatch, returncode=17)
    output = attempt_inputs["root"] / "engine-failure"

    with pytest.raises(RuntimeError, match="Sander exited 17"):
        _run(attempt_inputs, output)

    process = json.loads((output / "process.json").read_text())
    acceptance = json.loads((output / "acceptance.json").read_text())
    assert process["exit_status"] == 17
    assert process["ended_utc"] is not None
    assert acceptance["stage"] == "engine"
    assert "Sander exited 17" in acceptance["error"]


def test_ti1_writer_audits_qm_region_and_rejects_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "writer"
    output.mkdir()
    (output / "ti1.prmtop").write_text("topology\n")
    audited: list[Path] = []
    monkeypatch.setattr(
        qmmm,
        "audit_qm_region",
        lambda path: audited.append(path) or {"qmmask": qmmm.QM_MASK},
    )

    mdin = qmmm.write_ti1_dftb3_input(output, max_cycles=125)

    text = mdin.read_text()
    assert audited == [output / "ti1.prmtop"]
    assert f"qmmask='{qmmm.QM_MASK}'" in text
    assert "maxcyc=125" in text
    with pytest.raises(FileExistsError):
        qmmm.write_ti1_dftb3_input(output)


def test_legacy_ti1_runner_rejects_missing_explicit_arguments(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="requires explicit coordinates"):
        qmmm.run_ti1_dftb3(tmp_path, tmp_path)


def test_legacy_ti1_runner_forwards_explicit_arguments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    coordinates = tmp_path / "start.rst7"
    reference = tmp_path / "reference.rst7"
    checkpoint = tmp_path / "prior"
    output = tmp_path / "attempt"
    received: dict[str, object] = {}

    def fake_run_attempt(**kwargs: object) -> Path:
        received.update(kwargs)
        return output / "ti1_dftb3.mdout"

    monkeypatch.setattr(ti1_attempt, "run_attempt", fake_run_attempt)

    result = qmmm.run_ti1_dftb3(
        source,
        tmp_path / "amber" / "bin",
        max_cycles=125,
        coordinates=coordinates,
        reference=reference,
        attempt=output,
        target="poly_xama_maxa",
        checkpoint_attempt=checkpoint,
    )

    assert result == output / "ti1_dftb3.mdout"
    assert received == {
        "topology": source / "ti1.prmtop",
        "coordinates": coordinates,
        "reference": reference,
        "output": output,
        "amber_bin": tmp_path / "amber" / "bin",
        "target": "poly_xama_maxa",
        "max_cycles": 125,
        "checkpoint_attempt": checkpoint,
    }


def test_build_rejects_qmmm_before_creating_output(tmp_path: Path) -> None:
    output = tmp_path / "should-not-be-created"

    with pytest.raises(ValueError):
        build_nyl12_ti1_qmmm.build(
            target="pa6_l2",
            output=output,
            source_dir=tmp_path / "source",
            amber_bin=tmp_path / "amber" / "bin",
            run_qmmm=True,
        )

    assert not output.exists()
