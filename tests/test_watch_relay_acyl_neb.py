"""Check coherent NEB publication, terminal capture and exact rsync filter behavior."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pytest
from test_neb_progress import native_report

from nyl12_covalent import relay_neb_progress as helper
from nyl12_covalent.neb_progress import parse_last_complete_neb_report
from scripts import watch_relay_acyl_neb as monitor


def setup_monitor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, list]:
    """Build a24-image native report; fake only expensive reconstructed movie output."""
    case = tmp_path / "case.json"
    case.write_text(
        json.dumps(
            {
                "topology": "fixed.parm7",
                "qm_charge": -1,
                "qm_atom_ids": list(range(1, 103)),
                "amber_home": "/amber",
            }
        )
    )
    campaign = tmp_path / "04_neb"
    native = campaign / "calculation"
    native.mkdir(parents=True)
    config = campaign / "config.json"
    config.write_text(
        json.dumps(
            {
                "case_path": str(case),
                "topology": "fixed.parm7",
                "reactant": "reference.ncrst",
                "active_mask": "@1-4662",
                "active_atoms": 4662,
                "regular_images": 24,
            }
        )
    )
    (native / "minimized.mdout").write_bytes(
        native_report([-10 + i / 100 for i in range(24)])
        + b"Testing convergence in cycle 5\n"
    )
    (native / "nebpath.xyz").write_text("Synthetic fixture, reconstruction mocked\n")
    monkeypatch.setattr(monitor, "source_paths", lambda: [Path(__file__)])
    monkeypatch.setattr(monitor, "writer_state", lambda *_: (77, False))
    captures = []

    def capture(
        mdout: Path, xyz: Path, destination: Path, cfg: dict, count: int, **kwargs: Any
    ) -> Any:
        captures.append(kwargs)
        destination.mkdir()
        return parse_last_complete_neb_report(mdout.read_bytes(), count)

    def publish(
        cfg: dict,
        band: Path,
        report: Any,
        output: Path,
        identity: dict,
        terminal: bool,
        milestone: Path,
    ) -> dict:
        record = {
            "identity": identity,
            "terminal": terminal,
            "cycle": report.cycle,
            "latest_update": f"update_{len(captures):06d}",
        }
        helper.atomic_json(output / "progress.json", record)
        return record

    monkeypatch.setattr(monitor, "snapshot_complete_band_from_xyz", capture)
    monkeypatch.setattr(monitor, "publish_band", publish)
    return config, campaign / "progress", captures


def test_refresh_deduplicates_and_final_capture_requires_terminal_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg, output, captures = setup_monitor(tmp_path, monkeypatch)
    assert (
        monitor.refresh(cfg, output, tmp_path / "journal")["outcome"]
        == "published coherent band"
    )
    assert captures == [{"writer_pid": 77, "completed_job_status": None}]
    assert "unchanged" in monitor.refresh(cfg, output, tmp_path / "journal")["outcome"]
    assert len(captures) == 1
    assert len(list((output / "source_milestones").glob("milestone_*"))) == 1
    monkeypatch.setattr(monitor, "writer_state", lambda *_: (None, True))
    assert monitor.refresh(cfg, output, tmp_path / "journal")["terminal"]
    assert captures[-1]["completed_job_status"] == cfg.parent / "job/status.json"
    assert captures[-1]["writer_pid"] is None


def test_waiting_transient_and_fixed_input_failure_are_distinct(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg, output, _ = setup_monitor(tmp_path, monkeypatch)
    monkeypatch.setattr(monitor, "writer_state", lambda *_: (None, False))
    assert "waiting" in monitor.refresh(cfg, output, tmp_path / "journal")["outcome"]
    monkeypatch.setattr(monitor, "writer_state", lambda *_: (77, False))

    def moving(*args: Any, **kwargs: Any) -> None:
        raise RuntimeError("XYZ snapshot inputs changed during byte capture")

    monkeypatch.setattr(monitor, "snapshot_complete_band_from_xyz", moving)
    assert "changed" in monitor.refresh(cfg, output, tmp_path / "journal")["outcome"]

    def invalid(*args: Any, **kwargs: Any) -> None:
        raise RuntimeError("Fixed reference or topology changed while reconstructing")

    monkeypatch.setattr(monitor, "snapshot_complete_band_from_xyz", invalid)
    with pytest.raises(RuntimeError, match="Fixed reference"):
        monitor.refresh(cfg, output, tmp_path / "journal")
    config = json.loads(cfg.read_text())
    config["regular_images"] = 23
    cfg.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="24 regular"):
        monitor.refresh(cfg, output, tmp_path / "journal")


def test_parser_waits_for_first_complete_report_and_rejects_terminal_empty(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg, output, _ = setup_monitor(tmp_path, monkeypatch)
    mdout = cfg.parent / "calculation/minimized.mdout"
    mdout.write_bytes(b"NEB Report\nImg ")
    assert (
        "No complete" in monitor.refresh(cfg, output, tmp_path / "journal")["outcome"]
    )
    monkeypatch.setattr(monitor, "writer_state", lambda *_: (None, True))
    with pytest.raises(ValueError, match="No complete"):
        monitor.refresh(cfg, output, tmp_path / "journal")
    mdout.unlink()
    with pytest.raises(FileNotFoundError, match="Terminal NEB"):
        monitor.refresh(cfg, output, tmp_path / "journal")


def test_writer_state_uses_job_end_and_rejects_live_terminal_conflict(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(monitor, "find_sander_pid", lambda _: None)
    assert monitor.writer_state(tmp_path, Path("/sander")) == (None, False)
    status = tmp_path / "job/status.json"
    status.parent.mkdir()
    status.write_text(json.dumps({"ended_utc": "done", "exit_status": 0}))
    assert monitor.writer_state(tmp_path, Path("/sander")) == (None, True)
    monkeypatch.setattr(monitor, "find_sander_pid", lambda _: os.getpid())
    actual = Path(f"/proc/{os.getpid()}/exe").resolve()
    with pytest.raises(ValueError, match="executable differs"):
        monitor.writer_state(tmp_path, Path("/bin/bash"))
    with pytest.raises(RuntimeError, match="conflicts"):
        monitor.writer_state(tmp_path, actual)


def test_cli_watch_stops_after_terminal_and_records_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The CLI must poll at30s, stop after final export and durably report errors."""
    config = tmp_path / "config.json"
    config.write_text("{}")
    events = []
    sleeps = []
    monkeypatch.setattr(
        monitor,
        "audit",
        lambda journal, phase, action, files=None: events.append((phase, action)),
    )
    results = iter(
        [
            {"terminal": False, "outcome": "waiting"},
            {"terminal": True, "outcome": "published coherent band"},
        ]
    )
    monkeypatch.setattr(monitor, "refresh", lambda *_: next(results))
    monkeypatch.setattr(monitor.time, "sleep", lambda interval: sleeps.append(interval))
    monitor.main(["--config", str(config), "--watch", "--hours", "72"])
    assert sleeps == [30] and [e[0] for e in events] == [
        "intent",
        "outcome",
        "intent",
        "outcome",
    ]

    def fail(*_: Any) -> None:
        raise ValueError("invalid native band")

    failures = []
    monkeypatch.setattr(monitor, "refresh", fail)
    monkeypatch.setattr(
        monitor, "record_terminal_failure", lambda *args: failures.append(args)
    )
    with pytest.raises(ValueError, match="invalid native band"):
        monitor.main(["--config", str(config)])
    assert len(failures) == 1 and events[-1][0] == "outcome"
    with pytest.raises(SystemExit):
        monitor.main(["--config", str(config), "--hours", "nan"])
