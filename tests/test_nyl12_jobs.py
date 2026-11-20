"""Tests for the reproducible tmux job launcher."""

import json
import os
import subprocess
from pathlib import Path

import pytest

from nyl12_covalent.jobs import launch_job


def _fake_tmux(tmp_path: Path, occupied: bool = False) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir(parents=True)
    executable = bindir / "tmux"
    executable.write_text(
        "#!/bin/sh\n"
        "if [ \"$1\" = has-session ]; then exit "
        + ("0" if occupied else "1")
        + "; fi\nexit 0\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    return bindir


def test_launch_records_argv_and_hash_and_quotes_shell_metacharacters(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bindir = _fake_tmux(tmp_path)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    cwd = tmp_path / "work"
    cwd.mkdir()
    executable = cwd / "command"
    executable.write_text("#!/bin/sh\nprintf '%s' \"$1\"\n", encoding="utf-8")
    executable.chmod(0o755)
    input_path = cwd / "input.dat"
    input_path.write_text("contents", encoding="utf-8")
    output = tmp_path / "run"
    argv = [str(executable), "$(touch PWNED); spaces 'quote'"]
    script = launch_job("safe_session", output, cwd, argv, {"VALUE": "a b;$HOME"}, [input_path])

    recorded = json.loads((output / "command.json").read_text(encoding="utf-8"))
    assert recorded["command"] == argv
    assert recorded["input_sha256"][str(input_path.resolve())]
    process = subprocess.run(["bash", str(script)], cwd=cwd, check=False)
    assert process.returncode == 0
    assert not (cwd / "PWNED").exists()
    assert "$(touch PWNED); spaces 'quote'" in (output / "job.log").read_text()
    status = json.loads((output / "status.json").read_text(encoding="utf-8"))
    assert status["started_utc"].endswith("Z")
    assert status["ended_utc"].endswith("Z")
    assert status["exit_status"] == 0


def test_failure_status_is_retained(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bindir = _fake_tmux(tmp_path)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    cwd = tmp_path / "work"
    cwd.mkdir()
    output = tmp_path / "failed"
    script = launch_job("failed_session", output, cwd, ["/bin/sh", "-c", "echo failed; exit 7"])
    assert subprocess.run(["bash", str(script)], check=False).returncode == 7
    status = json.loads((output / "status.json").read_text(encoding="utf-8"))
    assert status["exit_status"] == 7
    assert "failed" in (output / "job.log").read_text(encoding="utf-8")


def test_rejects_existing_output_and_tmux_session(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bindir = _fake_tmux(tmp_path, occupied=True)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    cwd = tmp_path / "work"
    cwd.mkdir()
    with pytest.raises(ValueError, match="session already exists"):
        launch_job("busy", tmp_path / "new", cwd, ["/bin/true"])

    monkeypatch.setenv("PATH", f"{_fake_tmux(tmp_path / 'unused')}:{os.environ['PATH']}")
    output = tmp_path / "existing"
    output.mkdir()
    with pytest.raises(ValueError, match="output directory already exists"):
        launch_job("free", output, cwd, ["/bin/true"])


def test_rerun_uses_recorded_cwd_and_path_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bindir = _fake_tmux(tmp_path)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    cwd = tmp_path / "recorded workdir"
    cwd.mkdir()
    command_dir = tmp_path / "job tools"
    command_dir.mkdir()
    executable = command_dir / "report-cwd"
    executable.write_text("#!/bin/sh\npwd\n", encoding="utf-8")
    executable.chmod(0o755)
    output = tmp_path / "cwd-run"
    script = launch_job(
        "cwd_test", output, cwd, ["report-cwd"], {"PATH": str(command_dir)}
    )

    process = subprocess.run(
        ["bash", str(script)], cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert process.returncode == 0
    assert str(cwd) in (output / "job.log").read_text(encoding="utf-8")
    recorded = json.loads((output / "command.json").read_text(encoding="utf-8"))
    assert recorded["command"][0] == str(executable.resolve())


def test_rejects_tmux_session_name_injection(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bindir = _fake_tmux(tmp_path)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    cwd = tmp_path / "work"
    cwd.mkdir()
    with pytest.raises(ValueError, match="only letters"):
        launch_job("name:0", tmp_path / "output", cwd, ["/bin/true"])


def test_preserves_executable_symlink_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bindir = _fake_tmux(tmp_path)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    cwd = tmp_path / "work"
    cwd.mkdir()
    target = tmp_path / "real-python"
    target.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    target.chmod(0o755)
    link = tmp_path / "venv-python"
    link.symlink_to(target)
    output = tmp_path / "symlink-run"
    launch_job("symlink_test", output, cwd, [str(link)])
    recorded = json.loads((output / "command.json").read_text(encoding="utf-8"))
    assert recorded["command"][0] == str(link)


@pytest.mark.parametrize("explicit_override", [False, True])
def test_tmux_uses_system_libraries_and_job_retains_scientific_libraries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, explicit_override: bool
) -> None:
    """Keep Amber/Conda libtinfo away from tmux without losing the job setting."""
    bindir = _fake_tmux(tmp_path)
    executable = bindir / "tmux"
    executable.write_text(
        '#!/bin/sh\n'
        'if [ -n "${LD_LIBRARY_PATH-}" ]; then exit 86; fi\n'
        'if [ "$1" = has-session ]; then exit 1; fi\n'
        'exit 0\n'
    )
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    inherited = str(tmp_path / "parent_amber_lib")
    requested = str(tmp_path / "explicit_amber_lib")
    monkeypatch.setenv("LD_LIBRARY_PATH", inherited)
    cwd = tmp_path / "work"
    cwd.mkdir()
    output = tmp_path / "run"
    overrides = {"LD_LIBRARY_PATH": requested} if explicit_override else {}
    script = launch_job(
        "library_isolation", output, cwd,
        ["/bin/sh", "-c", 'printf "%s" "$LD_LIBRARY_PATH"'], overrides,
    )
    expected = requested if explicit_override else inherited
    recorded = json.loads((output / "command.json").read_text())
    assert recorded["env_overrides"]["LD_LIBRARY_PATH"] == expected
    subprocess.run(["bash", str(script)], check=True)
    assert (output / "job.log").read_text() == expected
