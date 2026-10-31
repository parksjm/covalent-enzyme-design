"""Tests for the portable Amber GPU queue runner."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_nyl12_md_queue.py"
SPEC = importlib.util.spec_from_file_location("run_nyl12_md_queue", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_run_queue_sets_gpu_and_runs_entries(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    campaign = tmp_path / "campaign"
    replica = campaign / "systems/wt/replica_01"
    replica.mkdir(parents=True)
    (replica / "run_replica.sh").write_text("exit 0\n")
    queue = tmp_path / "queue.txt"
    queue.write_text("systems/wt/replica_01\n")
    calls = []
    monkeypatch.setattr(
        MODULE.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs))
    )
    MODULE.run_queue(campaign, queue, 1)
    assert calls[0][1]["env"]["CUDA_VISIBLE_DEVICES"] == "1"
    assert calls[0][1]["cwd"] == replica


def test_run_queue_rejects_empty_escape_and_bad_gpu(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    queue = tmp_path / "queue.txt"
    queue.write_text("")
    with pytest.raises(ValueError, match="empty"):
        MODULE.run_queue(campaign, queue, 0)
    with pytest.raises(ValueError, match="nonnegative"):
        MODULE.run_queue(campaign, queue, -1)
    queue.write_text("../outside\n")
    with pytest.raises(ValueError, match="escapes"):
        MODULE.run_queue(campaign, queue, 0)
