"""Check preflight bookkeeping using a fake executable, never Gaussian."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "gaussian_preflight",
    Path(__file__).resolve().parents[1] / "scripts/preflight_gaussian16.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _installation(tmp_path: Path, message: str, status: int) -> tuple[Path, Path, Path]:
    install = tmp_path / "installation"
    profile = install / "g16/bsd/g16.profile"
    profile.parent.mkdir(parents=True)
    profile.write_text("# fake profile for unit test\n")
    executable = install / "g16/g16"
    executable.write_text(
        f"#!/bin/sh\ncat >/dev/null\necho '{message}'\nexit {status}\n"
    )
    executable.chmod(0o700)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    return install, scratch, tmp_path / "attempt"


def test_success_records_inputs_and_preserves_scratch(tmp_path: Path) -> None:
    install, scratch, output = _installation(
        tmp_path, "Normal termination of Gaussian 16", 0
    )
    MODULE.run_preflight(install, scratch, output)
    report = json.loads((output / "preflight.json").read_text())
    assert report["normal_termination"] is True
    assert report["exit_status"] == 0
    assert Path(report["scratch"]).is_dir()
    assert len(report["input_sha256"]) == 5
    assert "not a parameter calculation" in (output / "water.gjf").read_text()


@pytest.mark.parametrize(
    ("message", "status"),
    [("not finished", 0), ("Normal termination of Gaussian 16", 2)],
)
def test_failure_preserves_evidence(
    tmp_path: Path, message: str, status: int
) -> None:
    install, scratch, output = _installation(tmp_path, message, status)
    with pytest.raises(RuntimeError, match="preflight failed"):
        MODULE.run_preflight(install, scratch, output)
    report = json.loads((output / "preflight.json").read_text())
    assert report["exit_status"] == status
    assert (output / "water.log").is_file()


def test_existing_output_is_not_overwritten(tmp_path: Path) -> None:
    install, scratch, output = _installation(tmp_path, "unused", 0)
    output.mkdir()
    sentinel = output / "sentinel"
    sentinel.write_text("retain")
    with pytest.raises(FileExistsError):
        MODULE.run_preflight(install, scratch, output)
    assert sentinel.read_text() == "retain"


def test_inaccessible_executable_blocks_launch(tmp_path: Path) -> None:
    install, scratch, output = _installation(tmp_path, "unused", 0)
    (install / "g16/g16").chmod(0o600)
    with pytest.raises(PermissionError, match="inaccessible"):
        MODULE.run_preflight(install, scratch, output)
    assert not output.exists()
