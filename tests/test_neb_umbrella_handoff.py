"""Focused tests for procfs identity parsing and safe process handoff."""

import ctypes
import errno
import os
import signal

import pytest

from nyl12_covalent import neb_umbrella_handoff as handoff


def _stat_line(
    pid: int = 314,
    comm: str = "python (driver)",
    state: str = "S",
    ppid: int = 200,
    start_ticks: int = 12345,
    raw_exit_status: int | None = None,
) -> str:
    fields = ["0"] * 49
    fields[0] = state
    fields[1] = str(ppid)
    fields[19] = str(start_ticks)
    if raw_exit_status is not None:
        fields.extend(["0"] * (50 - len(fields)))
        fields[49] = str(raw_exit_status)
    return f"{pid} ({comm}) " + " ".join(fields)


def test_parse_stat_parenthesized_command_and_exit_statuses() -> None:
    for raw_status, expected in [(0, 0), (256, 1), (9, -9)]:
        parsed = handoff._parse_proc_stat(
            _stat_line(comm="sander (native)", state="Z", raw_exit_status=raw_status),
            314,
        )
        assert parsed["state"] == "Z"
        assert parsed["ppid"] == 200
        assert parsed["start_ticks"] == 12345
        assert os.waitstatus_to_exitcode(parsed["raw_exit_status"]) == expected


def test_parse_stat_rejects_pid_mismatch_and_missing_start_ticks() -> None:
    with pytest.raises(ValueError, match="PID"):
        handoff._parse_proc_stat(_stat_line(pid=315), 314)
    with pytest.raises(ValueError, match="start-time"):
        handoff._parse_proc_stat("314 (short) S 1", 314)


def _identity(pid: int = 314, **updates: object) -> dict[str, object]:
    result: dict[str, object] = {
        "pid": pid,
        "ppid": 200,
        "state": "S",
        "start_ticks": 12345,
        "argv": ["/usr/bin/python", "run_neb_umbrella.py", "window"],
        "cwd": "/campaign",
        "exe": "/usr/bin/python",
        "raw_exit_status": None,
    }
    result.update(updates)
    return result


def test_capture_reads_procfs_and_decodes_argv(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakePath:
        def __init__(self, value: str) -> None:
            self.value = value

        def __truediv__(self, item: str) -> "FakePath":
            return FakePath(f"{self.value}/{item}")

        def read_text(self, encoding: str) -> str:
            assert self.value.endswith("/stat") and encoding == "ascii"
            return _stat_line()

        def read_bytes(self) -> bytes:
            return b"python\0script.py\0\xff\0"

        def __str__(self) -> str:
            return self.value

    monkeypatch.setattr(handoff, "_proc_dir", lambda pid: FakePath(f"/proc/{pid}"))
    monkeypatch.setattr(handoff.os, "readlink", lambda p: f"{p}/target")
    captured = handoff.capture(314)
    assert captured["argv"] == ["python", "script.py", os.fsdecode(b"\xff")]
    assert captured["cwd"] == "/proc/314/cwd/target"
    assert captured["exe"] == "/proc/314/exe/target"


def test_identity_mismatch_rejects_pid_reuse(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(handoff, "capture", lambda pid: _identity(start_ticks=99999))
    with pytest.raises(RuntimeError, match="identity mismatch"):
        handoff._verify_identity(_identity())


def test_capture_missing_process_fails_explicitly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class MissingPath:
        def __truediv__(self, item: str) -> "MissingPath":
            return self

        def read_text(self, encoding: str) -> str:
            raise FileNotFoundError("missing stat")

    monkeypatch.setattr(handoff, "_proc_dir", lambda pid: MissingPath())
    with pytest.raises(ProcessLookupError):
        handoff.capture(314)


@pytest.mark.parametrize(
    ("native", "driver", "message"),
    [
        (
            _identity(exe="/usr/bin/sander", ppid=999, cwd="/campaign/stage"),
            _identity(pid=100),
            "parent",
        ),
        (
            _identity(exe="/usr/bin/sander", ppid=100, cwd="/campaign/stage"),
            _identity(pid=100, argv=["python", "script.py"]),
            "controller",
        ),
        (
            _identity(exe="/usr/bin/other", ppid=100, cwd="/campaign/stage"),
            _identity(pid=100),
            "sander",
        ),
        (
            _identity(exe="/usr/bin/sander", ppid=100, cwd="/elsewhere/stage"),
            _identity(pid=100),
            "campaign window",
        ),
    ],
)
def test_suspend_driver_rejects_invalid_relationship(
    monkeypatch: pytest.MonkeyPatch,
    native: dict[str, object],
    driver: dict[str, object],
    message: str,
) -> None:
    monkeypatch.setattr(
        handoff, "capture", lambda pid: native if pid == 314 else driver
    )
    monkeypatch.setattr(
        handoff, "signal_verified", lambda *_: pytest.fail("must not signal")
    )
    with pytest.raises(ValueError, match=message):
        handoff.suspend_driver(native, driver)


def test_native_completion_handles_running_and_zombie(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identity = _identity(exe="/usr/bin/sander", argv=["sander"])
    current = dict(identity)
    monkeypatch.setattr(handoff, "capture", lambda pid: current)
    assert handoff.native_completion(identity) is None
    current.update(state="Z", raw_exit_status=256, argv=[], cwd=None, exe=None)
    assert handoff.native_completion(identity) == 1
    current["start_ticks"] = 999
    with pytest.raises(RuntimeError, match="zombie process identity mismatch"):
        handoff.native_completion(identity)


def test_signal_verified_uses_pidfd_and_rejects_reused_process(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identity = _identity()
    sent: list[tuple[int, int]] = []
    monkeypatch.setattr(handoff, "_open_pidfd", lambda pid: 88)
    monkeypatch.setattr(
        handoff, "_send_pidfd_signal", lambda fd, sig: sent.append((fd, sig))
    )
    monkeypatch.setattr(handoff, "capture", lambda pid: _identity(start_ticks=999))
    monkeypatch.setattr(handoff.os, "close", lambda fd: None)
    with pytest.raises(RuntimeError, match="identity mismatch"):
        handoff.signal_verified(identity, signal.SIGTERM)
    assert sent == []


def test_release_driver_sends_term_then_cont_on_one_pidfd(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identity = _identity()
    sent: list[tuple[int, int]] = []
    closed: list[int] = []
    monkeypatch.setattr(handoff, "_open_pidfd", lambda pid: 88)
    monkeypatch.setattr(
        handoff, "_send_pidfd_signal", lambda fd, sig: sent.append((fd, sig))
    )
    monkeypatch.setattr(handoff, "capture", lambda pid: identity)
    monkeypatch.setattr(handoff.os, "close", closed.append)
    monkeypatch.setattr(
        handoff.os,
        "killpg",
        lambda *_: pytest.fail("process-group signaling is forbidden"),
    )
    handoff.release_driver(identity)
    assert sent == [(88, signal.SIGTERM), (88, signal.SIGCONT)]
    assert closed == [88]


def test_suspend_driver_accepts_campaign_window_stage(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native = _identity(
        exe="/usr/bin/sander",
        argv=["sander", "-O"],
        ppid=100,
        cwd="/campaign/window_044/03_production_block_04_5ps",
    )
    driver = _identity(
        pid=100,
        argv=["python", "run_neb_umbrella.py", "window"],
        cwd="/campaign",
    )
    stopped_driver = dict(driver)
    stopped_driver["state"] = "T"
    sent: list[tuple[int, int]] = []
    monkeypatch.setattr(
        handoff,
        "capture",
        lambda pid: native if pid == native["pid"] else stopped_driver,
    )
    monkeypatch.setattr(
        handoff,
        "signal_verified",
        lambda identity, sig: sent.append((identity["pid"], sig)),
    )
    handoff.suspend_driver(native, driver)
    assert sent == [(100, signal.SIGSTOP)]


def test_capture_zombie_accepts_missing_command_and_link_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class ZombiePath:
        def __init__(self, value: str) -> None:
            self.value = value

        def __truediv__(self, item: str) -> "ZombiePath":
            return ZombiePath(f"{self.value}/{item}")

        def read_text(self, encoding: str) -> str:
            return _stat_line(state="Z", raw_exit_status=0)

        def read_bytes(self) -> bytes:
            raise FileNotFoundError("zombie cmdline")

        def __str__(self) -> str:
            return self.value

    monkeypatch.setattr(handoff, "_proc_dir", lambda pid: ZombiePath(f"/proc/{pid}"))

    def missing_link(path: object) -> str:
        raise FileNotFoundError("zombie link")

    monkeypatch.setattr(handoff.os, "readlink", missing_link)
    result = handoff.capture(314)
    assert result["state"] == "Z"
    assert result["argv"] == []
    assert result["cwd"] is None
    assert result["exe"] is None
    assert result["raw_exit_status"] == 0


def test_libc_pidfd_bindings_send_and_report_system_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class FakeFunction:
        def __init__(self, result: int) -> None:
            self.result = result
            self.calls: list[tuple[object, ...]] = []
            self.argtypes: object = None
            self.restype: object = None

        def __call__(self, *args: object) -> int:
            self.calls.append(args)
            return self.result

    open_function = FakeFunction(73)
    send_function = FakeFunction(0)

    class FakeLibc:
        pidfd_open = open_function
        pidfd_send_signal = send_function

    monkeypatch.setattr(handoff, "_LIBC", FakeLibc())
    assert handoff._open_pidfd(512) == 73
    handoff._send_pidfd_signal(73, signal.SIGTERM)
    assert open_function.calls == [(512, 0)]
    assert send_function.calls == [(73, signal.SIGTERM, None, 0)]

    def denied_open(pid: int, flags: int) -> int:
        ctypes.set_errno(errno.ESRCH)
        return -1

    class DeniedLibc:
        pidfd_open = staticmethod(denied_open)

    monkeypatch.setattr(handoff, "_LIBC", DeniedLibc())
    with pytest.raises(ProcessLookupError):
        handoff._open_pidfd(512)

    def denied_send(pidfd: int, sig: int, info: object, flags: int) -> int:
        ctypes.set_errno(errno.EPERM)
        return -1

    class DeniedSignalLibc:
        pidfd_send_signal = staticmethod(denied_send)

    monkeypatch.setattr(handoff, "_LIBC", DeniedSignalLibc())
    with pytest.raises(PermissionError):
        handoff._send_pidfd_signal(73, signal.SIGTERM)


def test_pidfd_bindings_fail_when_libc_symbols_are_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(handoff, "_LIBC", object())
    with pytest.raises(RuntimeError, match="libc pidfd_open"):
        handoff._open_pidfd(512)
    with pytest.raises(RuntimeError, match="libc pidfd_send_signal"):
        handoff._send_pidfd_signal(73, signal.SIGTERM)
