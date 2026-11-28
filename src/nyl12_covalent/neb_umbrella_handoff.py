"""Safely hand off stopped umbrella-window drivers after native Sander blocks."""

from __future__ import annotations

import ctypes
import os
import signal
import time
from pathlib import Path
from typing import Any


_IDENTITY_FIELDS = ("start_ticks", "argv", "cwd", "exe")
_LIBC = ctypes.CDLL(None, use_errno=True)


def _parse_proc_stat(contents: str, expected_pid: int) -> dict[str, Any]:
    """Parse selected fields from Linux proc stat, including parenthesized comm."""
    opening = contents.find(" (")
    closing = contents.rfind(")")
    if opening <= 0 or closing <= opening:
        raise ValueError("malformed /proc stat command field")
    pid_text = contents[:opening]
    if not pid_text.isdecimal() or int(pid_text) != expected_pid:
        raise ValueError("/proc stat PID does not match requested PID")
    fields = contents[closing + 1 :].split()
    if len(fields) < 20:
        raise ValueError("/proc stat lacks required start-time field")
    try:
        ppid = int(fields[1])
        start_ticks = int(fields[19])
        raw_exit_status = int(fields[49]) if len(fields) >= 50 else None
    except ValueError as exc:
        raise ValueError("/proc stat contains an invalid numeric field") from exc
    if len(fields[0]) != 1:
        raise ValueError("/proc stat contains an invalid process state")
    return {
        "pid": expected_pid,
        "ppid": ppid,
        "state": fields[0],
        "start_ticks": start_ticks,
        "raw_exit_status": raw_exit_status,
    }


def _proc_dir(pid: int) -> Path:
    return Path("/proc") / str(pid)


def capture(pid: int) -> dict[str, Any]:
    """Capture a Linux process identity and its current state.

    Args:
        pid: Linux process ID to inspect.

    Returns:
        Process identity fields plus the raw wait status when procfs exposes it.
        Linux may omit argv, cwd and exe after a process becomes a zombie; those
        fields are then empty or ``None`` while stat identity remains available.

    Raises:
        ProcessLookupError: If the process disappears during capture.
        OSError: If procfs metadata cannot be read.
        ValueError: If procfs returns malformed data.
    """
    directory = _proc_dir(pid)
    try:
        stat = _parse_proc_stat((directory / "stat").read_text(encoding="ascii"), pid)
    except FileNotFoundError as exc:
        raise ProcessLookupError(pid, "process disappeared during capture") from exc

    try:
        argv_data = (directory / "cmdline").read_bytes()
    except FileNotFoundError as exc:
        if stat["state"] != "Z":
            raise ProcessLookupError(pid, "process disappeared during capture") from exc
        argv_data = b""
    argv = [os.fsdecode(item) for item in argv_data.split(b"\0") if item]

    def read_link(name: str) -> str | None:
        try:
            return os.readlink(directory / name)
        except FileNotFoundError as exc:
            if stat["state"] != "Z":
                raise ProcessLookupError(
                    pid, "process disappeared during capture"
                ) from exc
            return None

    cwd = read_link("cwd")
    exe = read_link("exe")
    return {**stat, "argv": argv, "cwd": cwd, "exe": exe}


def _require_identity(identity: dict[str, Any]) -> dict[str, Any]:
    try:
        pid = identity["pid"]
        if not isinstance(pid, int) or pid <= 0:
            raise ValueError
        for field in _IDENTITY_FIELDS:
            identity[field]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("incomplete process identity") from exc
    return identity


def _verify_identity(identity: dict[str, Any]) -> None:
    expected = _require_identity(identity)
    current = capture(expected["pid"])
    if any(current[field] != expected[field] for field in _IDENTITY_FIELDS):
        raise RuntimeError(f"process identity mismatch for PID {expected['pid']}")


def _open_pidfd(pid: int) -> int:
    """Open a pidfd through libc, failing explicitly when unavailable."""
    try:
        function = _LIBC.pidfd_open
    except AttributeError as exc:
        raise RuntimeError("libc pidfd_open is required") from exc
    function.argtypes = (ctypes.c_int, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    pidfd = function(pid, 0)
    if pidfd < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    return pidfd


def _send_pidfd_signal(pidfd: int, sig: int) -> None:
    """Send one signal through libc's pidfd_send_signal binding."""
    try:
        function = _LIBC.pidfd_send_signal
    except AttributeError as exc:
        raise RuntimeError("libc pidfd_send_signal is required") from exc
    function.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    result = function(pidfd, sig, None, 0)
    if result < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))


def _pidfd_send(identity: dict[str, Any], signals: tuple[int, ...]) -> None:
    """Verify once against a pidfd and send signals only to that process."""
    expected = _require_identity(identity)
    pidfd = _open_pidfd(expected["pid"])
    try:
        _verify_identity(expected)
        for sig in signals:
            _send_pidfd_signal(pidfd, sig)
    finally:
        os.close(pidfd)


def signal_verified(identity: dict[str, Any], sig: int) -> None:
    """Send a signal to a process only after verifying its pidfd identity.

    Args:
        identity: Previously captured process identity.
        sig: Signal number to send to that process.

    Raises:
        RuntimeError: If pidfd support is unavailable or identity changed.
        OSError: If opening the pidfd or sending the signal fails.
    """
    _pidfd_send(identity, (sig,))


def native_completion(identity: dict[str, Any]) -> int | None:
    """Return a child's exit code after it becomes a zombie, or None if running.

    Args:
        identity: Previously captured native process identity.

    Returns:
        Normalized exit code for a zombie process; otherwise None. For zombies,
        Linux may have removed argv/cwd/exe, so PID, start ticks and parent PID
        are checked against the earlier live identity.

    Raises:
        RuntimeError: If identity changed or zombie exit status is unavailable.
        ProcessLookupError: If the native process disappeared before inspection.
    """
    expected = _require_identity(identity)
    current = capture(expected["pid"])
    if current["state"] != "Z":
        if any(current[field] != expected[field] for field in _IDENTITY_FIELDS):
            raise RuntimeError(f"process identity mismatch for PID {expected['pid']}")
        return None
    if (
        current["start_ticks"] != expected["start_ticks"]
        or current["ppid"] != expected["ppid"]
    ):
        raise RuntimeError(
            f"zombie process identity mismatch for PID {expected['pid']}"
        )
    status = current["raw_exit_status"]
    if status is None:
        raise RuntimeError("zombie process has no procfs exit status")
    return os.waitstatus_to_exitcode(status)


def _is_sander(native: dict[str, Any]) -> bool:
    return Path(native["exe"]).name == "sander" and bool(native["argv"])


def suspend_driver(native: dict[str, Any], driver: dict[str, Any]) -> None:
    """Stop a validated umbrella window driver while its native block runs.

    The native Sander process is never signaled. The driver must be the native
    process's parent, run ``run_neb_umbrella.py`` with a literal ``window``
    argument, and have its campaign directory as the native stage's parent.

    Args:
        native: Captured native Sander identity.
        driver: Captured window-driver Python identity.

    Raises:
        ValueError: If the parent, command, or campaign/stage layout is invalid.
        RuntimeError: If either identity changes or the driver fails to stop.
        OSError: If process inspection or signaling fails.
    """
    _verify_identity(native)
    _verify_identity(driver)
    if driver["pid"] != native["ppid"]:
        raise ValueError("driver is not the native Sander parent")
    argv = driver["argv"]
    if "window" not in argv or not any(
        Path(argument).name == "run_neb_umbrella.py" for argument in argv
    ):
        raise ValueError("driver command is not an umbrella window controller")
    if not _is_sander(native):
        raise ValueError("native process executable is not sander")
    stage = Path(native["cwd"])
    window = stage.parent
    if (
        window.parent != Path(driver["cwd"])
        or not window.name.startswith("window_")
        or not window.name.removeprefix("window_").isdecimal()
    ):
        raise ValueError("native stage is not under a campaign window directory")

    signal_verified(driver, signal.SIGSTOP)
    deadline = time.monotonic() + 2.0
    while time.monotonic() < deadline:
        if capture(driver["pid"])["state"] in {"T", "t"}:
            return
        time.sleep(0.01)
    raise RuntimeError("window driver did not enter a stopped state")


def release_driver(driver: dict[str, Any]) -> None:
    """Terminate then continue a stopped driver after native exit was observed.

    The caller must first observe native Sander completion with
    :func:`native_completion`. SIGTERM becomes pending while the Python driver
    is stopped; SIGCONT lets it exit before it can launch another native block.
    Both signals use the same verified pidfd. The native process is never sent
    a signal.

    Args:
        driver: Captured identity of the stopped window-driver process.

    Raises:
        RuntimeError: If pidfd support is unavailable or identity changed.
        OSError: If process inspection or signaling fails.
    """
    _pidfd_send(driver, (signal.SIGTERM, signal.SIGCONT))
