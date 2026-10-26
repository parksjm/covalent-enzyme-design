"""Tests for the phase-II archive and resume utility."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.archive_invalid_phase2 import STAGE1, archive_phase2


def _make_run(root: Path) -> Path:
    """Create a valid synthetic run with phase-II artifacts."""
    root.mkdir()
    (root / "status.json").write_text(
        json.dumps({"completed_windows": STAGE1 + ["stage2_00", "stage2_06"]})
    )
    for name in STAGE1:
        (root / name).mkdir()
    stage = root / "stage1_20"
    (stage / "acceptance.json").write_text(
        json.dumps({"distances_A": {"OG1_HG1": 1.31, "Asp272_OD1_HG1": 1.19}})
    )
    for name in ("stage2_00", "stage2_06"):
        directory = root / name
        directory.mkdir()
        (directory / "result.dat").write_text(name)
    return root


def test_dry_run_does_not_mutate(tmp_path: Path) -> None:
    """Dry run validates inputs but leaves run and clutter untouched."""
    run = _make_run(tmp_path / "run")
    archive = archive_phase2(run, tmp_path / "clutter", dry_run=True)
    assert not archive.exists()
    assert (run / "stage2_06").exists()
    assert not (archive / "status_before.json").exists()


def test_archives_all_phase2_and_resets_status(tmp_path: Path) -> None:
    """Archive incomplete phase-II windows and preserve file hashes."""
    run = _make_run(tmp_path / "run")
    archive = archive_phase2(run, tmp_path / "clutter")
    assert {path.name for path in archive.glob("stage2_*")} == {
        "stage2_00",
        "stage2_06",
    }
    assert not (run / "stage2_06").exists()
    assert (archive / "status_before.json").exists()
    status = json.loads((run / "status.json").read_text())
    assert status["completed_windows"] == STAGE1
    assert status["state"] == "stopped"
    assert "donor-separation wall" in status["reason"]
    manifest = json.loads((archive / "manifest.json").read_text())
    assert all(
        record["sha256_before"] == record["sha256_after"]
        for record in manifest["directories"]
    )


@pytest.mark.parametrize(("og_hg", "asp_hg"), [(1.3, 1.0), (1.4, 1.2), (1.0, 1.0)])
def test_rejects_invalid_acceptance(
    tmp_path: Path, og_hg: float, asp_hg: float
) -> None:
    """Reject geometry that fails either strict acceptance threshold."""
    run = _make_run(tmp_path / "run")
    (run / "stage1_20" / "acceptance.json").write_text(
        json.dumps({"distances_A": {"OG1_HG1": og_hg, "Asp272_OD1_HG1": asp_hg}})
    )
    with pytest.raises(ValueError, match="geometry gates"):
        archive_phase2(run, tmp_path / "clutter")
    assert (run / "stage2_06").exists()
