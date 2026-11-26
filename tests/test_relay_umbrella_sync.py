"""Compact sync transfers allowlisted artifacts and holds readers on failure."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

import pytest
from scripts import sync_relay_umbrella as sync


def config(tmp_path: Path) -> argparse.Namespace:
    """Return one explicit temporary sync configuration."""
    return argparse.Namespace(
        campaign=tmp_path,
        execution=tmp_path,
        identity=tmp_path / "key",
        journal_root=tmp_path,
        benchmark_source=None,
        remote="azalea",
        remote_python="/python",
        inspector=Path("/inspect.py"),
    )


def test_compact_sync_allowlist_and_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = []
    monkeypatch.setattr(sync, "journal", lambda *args: None)
    monkeypatch.setattr(
        sync.subprocess, "run", lambda command, **kwargs: calls.append(command)
    )
    sync.sync_once(config(tmp_path))
    assert not (tmp_path / ".compact_sync_in_progress").exists()
    command = calls[-1]
    assert "umbrella_coordinate.dat" in command
    assert "*.csv" in command and "*.npz" in command
    assert "*.nc" not in command and "*.ncrst" not in command
    assert command.index("failed_analysis_*/***") < command.index("*/")
    assert command[-3] == "*"


def test_sync_failure_retains_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sync, "journal", lambda *args: None)

    def fail(*args: object, **kwargs: object) -> None:
        raise subprocess.CalledProcessError(1, ["ssh"])

    monkeypatch.setattr(sync.subprocess, "run", fail)
    with pytest.raises(subprocess.CalledProcessError):
        sync.sync_once(config(tmp_path))
    assert (tmp_path / ".compact_sync_in_progress").is_file()
