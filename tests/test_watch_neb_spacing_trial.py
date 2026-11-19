"""Verify isolated spacing monitors retain own source and honor24/48 image counts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from test_neb_progress import native_report

from nyl12_covalent.neb_progress import parse_last_complete_neb_report
from nyl12_covalent.relay_neb_progress import atomic_json
from scripts import watch_neb_spacing_trial as monitor


@pytest.mark.parametrize("count", [24, 48])
def test_spacing_watcher_passes_exact_band_count_and_retains_own_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, count: int
) -> None:
    """Parser/capture/export must all use the declared regular count, including48."""
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
    native = tmp_path / "calculation"
    native.mkdir()
    (native / "minimized.mdout").write_bytes(
        native_report([-10 + i / 100 for i in range(count)])
        + b"Testing convergence in cycle 5\n"
    )
    (native / "nebpath.xyz").write_text(
        "Native capture mocked after count validation\n"
    )
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "regular_images": count,
                "active_atoms": 4662,
                "case_path": str(case),
                "topology": "fixed.parm7",
            }
        )
    )
    output = tmp_path / "progress"
    calls = []
    monkeypatch.setattr(monitor, "source_paths", lambda: [Path(__file__)])
    monkeypatch.setattr(monitor, "writer_state", lambda *_: (77, False))

    def capture(
        mdout: Path,
        xyz: Path,
        destination: Path,
        cfg: dict,
        expected: int,
        **kwargs: Any,
    ) -> Any:
        assert expected == count == cfg["regular_images"]
        assert kwargs == {"writer_pid": 77, "completed_job_status": None}
        calls.append(expected)
        destination.mkdir()
        return parse_last_complete_neb_report(mdout.read_bytes(), expected)

    def publish(
        cfg: dict,
        band: Path,
        report: Any,
        folder: Path,
        identity: dict,
        terminal: bool,
        milestone: Path,
    ) -> dict:
        assert len(report.images) == count
        assert str(Path(monitor.__file__).resolve()) in identity["source_sha256"]
        assert any(
            p.name.endswith("watch_neb_spacing_trial.py") for p in milestone.iterdir()
        )
        record = {
            "identity": identity,
            "terminal": terminal,
            "cycle": report.cycle,
            "latest_update": "update_000001",
        }
        atomic_json(folder / "progress.json", record)
        return record

    monkeypatch.setattr(monitor, "snapshot_complete_band_from_xyz", capture)
    monkeypatch.setattr(monitor, "publish_band", publish)
    assert (
        monitor.refresh(config, output, tmp_path / "journal")["outcome"]
        == "published coherent band"
    )
    assert calls == [count]
    assert (
        "unchanged" in monitor.refresh(config, output, tmp_path / "journal")["outcome"]
    )
    assert calls == [count]


@pytest.mark.parametrize("count", [0, 25, 49])
def test_spacing_watcher_rejects_unsupported_regular_counts(
    tmp_path: Path, count: int
) -> None:
    """The native extra slot must not be confused with a supported regular count."""
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"regular_images": count, "active_atoms": 4662}))
    with pytest.raises(ValueError, match="24 or48"):
        monitor.refresh(config, tmp_path / "progress", tmp_path / "journal")
