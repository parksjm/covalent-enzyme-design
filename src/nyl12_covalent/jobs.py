"""Create reproducible, logged jobs in detached tmux sessions."""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Mapping, Sequence


def launch_job(
    session_name: str,
    output_dir: Path,
    cwd: Path,
    command: Sequence[str],
    env_overrides: Mapping[str, str] | None = None,
    input_paths: Sequence[Path] = (),
) -> Path:
    """Create a provenance bundle and launch a command in detached tmux.

    Args:
        session_name: Unique tmux session name.
        output_dir: New directory for scripts, metadata, and logs.
        cwd: Existing working directory for the command.
        command: Executable and argument vector, passed without shell parsing.
        env_overrides: Explicit environment values for the job.
        input_paths: Optional files whose SHA-256 hashes should be recorded.

    Returns:
        Path to the generated run script.

    Raises:
        ValueError: If inputs are invalid or the session/output already exists.
        FileNotFoundError: If the cwd, executable, input, or tmux is missing.
        RuntimeError: If tmux cannot start the session.
    """
    allowed_session_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
    if not session_name or not all(char in allowed_session_chars for char in session_name):
        raise ValueError("session_name may contain only letters, digits, '_' and '-'")
    if not command or not command[0]:
        raise ValueError("command must contain an executable")
    workdir = cwd.resolve(strict=True)
    if not workdir.is_dir():
        raise ValueError(f"cwd is not a directory: {workdir}")
    overrides = dict(env_overrides or {})
    for key, value in overrides.items():
        if not key.isidentifier() or "\x00" in value:
            raise ValueError(f"invalid environment override: {key!r}")
    executable = command[0]
    if os.sep in executable:
        candidate = Path(executable)
        resolved_executable = Path(os.path.abspath(workdir / candidate if not candidate.is_absolute() else candidate))
    else:
        found = shutil.which(executable, path=overrides.get("PATH", os.environ.get("PATH")))
        resolved_executable = Path(found) if found else None
    if (
        not resolved_executable
        or not resolved_executable.is_file()
        or not os.access(resolved_executable, os.X_OK)
    ):
        raise FileNotFoundError(f"executable not found or not executable: {executable}")
    tmux = shutil.which("tmux")
    if tmux is None:
        raise FileNotFoundError("tmux executable not found")
    destination = output_dir.expanduser().absolute()
    if destination.exists():
        raise ValueError(f"output directory already exists: {destination}")
    if subprocess.run([tmux, "has-session", "-t", f"={session_name}"], check=False).returncode == 0:
        raise ValueError(f"tmux session already exists: {session_name}")

    inputs = {str(path.resolve(strict=True)): _sha256(path) for path in input_paths}
    destination.mkdir(parents=True)
    log_path = destination / "job.log"
    status_path = destination / "status.json"
    run_path = destination / "run.sh"
    metadata = {
        "session_name": session_name,
        "cwd": str(workdir),
        "command": [str(resolved_executable), *command[1:]],
        "env_overrides": overrides,
        "input_sha256": inputs,
        "started_utc": None,
        "ended_utc": None,
        "exit_status": None,
        "log": str(log_path),
    }
    (destination / "command.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (destination / "status.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    run_path.write_text(
        _render_script(
            [str(resolved_executable), *command[1:]], overrides, workdir,
            log_path, status_path, Path(os.path.abspath(sys.executable))
        ),
        encoding="utf-8",
    )
    run_path.chmod(0o750)
    result = subprocess.run(
        [tmux, "new-session", "-d", "-s", session_name, "-c", str(workdir),
         shlex.join(["bash", str(run_path)])],
        check=False,
    )
    if result.returncode:
        raise RuntimeError(f"tmux failed to launch session {session_name!r}")
    return run_path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _render_script(
    command: Sequence[str],
    overrides: Mapping[str, str],
    cwd: Path,
    log_path: Path,
    status_path: Path,
    python: Path,
) -> str:
    assignments = "\n".join(f"export {key}={shlex.quote(value)}" for key, value in overrides.items())
    updater = (
        "import datetime,json,os; p=os.environ['STATUS_PATH']; d=json.load(open(p)); "
        "d['started_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('+00:00','Z'); "
        "json.dump(d,open(p,'w'),indent=2,sort_keys=True); open(p,'a').write('\\n')"
    )
    finalizer = (
        "import datetime,json,os,tempfile; p=os.environ['STATUS_PATH']; d=json.load(open(p)); "
        "d['ended_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat().replace('+00:00','Z'); "
        "d['exit_status']=int(os.environ['JOB_EXIT']); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(p)); "
        "f=os.fdopen(fd,'w'); json.dump(d,f,indent=2,sort_keys=True); f.write('\\n'); f.close(); os.replace(tmp,p)"
    )
    start_command = shlex.join([str(python), "-c", updater])
    finish_command = shlex.join([str(python), "-c", finalizer])
    return (
        "#!/usr/bin/env bash\n"
        f"export STATUS_PATH={shlex.quote(str(status_path))}\n"
        f"exec >>{shlex.quote(str(log_path))} 2>&1\n"
        f"cd -- {shlex.quote(str(cwd))} || exit $?\n"
        f"{assignments}\n"
        f"{start_command} || exit $?\n"
        "finish_job() { JOB_EXIT=$?; export JOB_EXIT; trap - EXIT HUP INT TERM; "
        f"{finish_command} || exit 125; exit \"$JOB_EXIT\"; }}\n"
        "trap finish_job EXIT\n"
        "trap 'exit 129' HUP\ntrap 'exit 130' INT\ntrap 'exit 143' TERM\n"
        f"{shlex.join(list(command))}\n"
    )
