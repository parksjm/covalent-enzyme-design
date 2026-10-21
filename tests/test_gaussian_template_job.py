"""Exercise Gaussian execution, failed acceptance, and cache lifecycle."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import run_nyl12_gaussian_template as worker


@pytest.fixture
def setup_job(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    source, scratch, install, amber, root = [
        tmp_path / name for name in ("prepared", "scratch", "install", "amber", "repo")
    ]
    for directory in (
        source,
        scratch,
        install / "g16/bsd",
        amber,
        root / "src/nyl12_covalent",
    ):
        directory.mkdir(parents=True)
    for name in ("input.gjf", "molecule.sdf", "atom_map.json"):
        (source / name).write_text(name)
    for path in (install / "g16/g16", amber / "espgen"):
        path.write_text("#!/bin/sh\nexit 0\n")
        path.chmod(0o750)
    (install / "g16/bsd/g16.profile").write_text("# synthetic profile\n")
    for name in ("gaussian_checks.py", "parameter_cache.py", "terminal_template.py"):
        (root / "src/nyl12_covalent" / name).write_text(name)
    (root / "uv.lock").write_text("fixture lock")
    monkeypatch.setattr(worker, "ROOT", root)
    monkeypatch.setattr(
        worker,
        "assess_gaussian",
        lambda *_: {"accepted": True, "final_coordinates_A": [[0, 0, 0]]},
    )
    return {
        "source": source,
        "scratch_root": scratch,
        "install_root": install,
        "amber_bin": amber,
        "cache_root": tmp_path / "cache",
        "output": tmp_path / "attempt",
    }


def _engine(monkeypatch: pytest.MonkeyPatch, code: int = 0) -> list[list[str]]:
    calls = []

    def execute(command: list[str], **kwargs: object) -> SimpleNamespace:
        calls.append(command)
        if command[0] == "bash":
            Path(kwargs["cwd"], "gaussian.log").write_text("synthetic output")
        else:
            Path(command[-1]).write_text("1 10\nESP data\n")
        return SimpleNamespace(returncode=code)

    monkeypatch.setattr(worker.subprocess, "run", execute)
    return calls


def test_success_cached_then_reused_without_engine(
    setup_job: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = _engine(monkeypatch)
    worker.run_template(**setup_job)
    assert len(calls) == 2
    acceptance = json.loads((setup_job["output"] / "acceptance.json").read_text())
    assert acceptance["accepted"] and not acceptance["fitted_forcefield_ready"]
    setup_job["output"] = setup_job["output"].with_name("reuse")
    worker.run_template(**setup_job)
    assert len(calls) == 2
    assert (setup_job["output"] / "cache_reuse.json").is_file()


@pytest.mark.parametrize("failure", ["engine", "science", "esp"])
def test_failure_keeps_records_and_does_not_cache(
    setup_job: dict[str, Path], monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    _engine(monkeypatch, 17 if failure == "engine" else 0)
    if failure == "science":
        monkeypatch.setattr(worker, "assess_gaussian", lambda *_: {"accepted": False})
    elif failure == "esp":
        monkeypatch.setattr(
            worker,
            "assess_gaussian",
            lambda *_: {"accepted": True, "final_coordinates_A": []},
        )
    with pytest.raises((RuntimeError, ValueError)):
        worker.run_template(**setup_job)
    assert not json.loads((setup_job["output"] / "acceptance.json").read_text())[
        "accepted"
    ]
    assert json.loads((setup_job["output"] / "process.json").read_text())["ended_utc"]
    assert not list(setup_job["cache_root"].glob("*/manifest.json"))


def test_rejects_existing_output(setup_job: dict[str, Path]) -> None:
    setup_job["output"].mkdir()
    with pytest.raises(FileExistsError):
        worker.run_template(**setup_job)


@pytest.mark.parametrize("missing", ["g16", "espgen", "scratch"])
def test_prerequisites_explicit(setup_job: dict[str, Path], missing: str) -> None:
    if missing == "g16":
        (setup_job["install_root"] / "g16/g16").chmod(0o600)
    elif missing == "espgen":
        (setup_job["amber_bin"] / "espgen").chmod(0o600)
    else:
        setup_job["scratch_root"] = setup_job["source"] / "input.gjf"
    with pytest.raises(PermissionError):
        worker.run_template(**setup_job)


def test_cli_forwards_explicit_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    captured = []
    monkeypatch.setattr(worker, "run_template", lambda *args: captured.extend(args))
    arguments = ["worker"]
    names = (
        "source",
        "output",
        "install-root",
        "scratch-root",
        "amber-bin",
        "cache-root",
    )
    for name in names:
        arguments.extend(["--" + name, "/" + name])
    monkeypatch.setattr("sys.argv", arguments)
    worker.main()
    assert captured == [Path("/" + name) for name in names]
