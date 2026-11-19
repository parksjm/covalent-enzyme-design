"""Verify the dependency-free NEB profile watcher and candidate provenance."""

import importlib.util
import json
import xml.etree.ElementTree as ET
from pathlib import Path

from test_neb_progress import native_report

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "watch_neb_progress.py"
_SPEC = importlib.util.spec_from_file_location("watch_neb_progress", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
watcher = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(watcher)


def test_svg_archive_and_candidates_are_coherent(tmp_path: Path) -> None:
    source = tmp_path / "neb.mdout"
    source.write_bytes(
        native_report([-10, -8, -9, -7, -11], (4, -6.5))
        + b"Testing convergence in cycle 3\n"
    )
    output = tmp_path / "watch"
    report = watcher.update_progress_once(source, output, 5)
    root = ET.fromstring((output / "neb_profile.svg").read_text())
    assert root.tag.endswith("svg")
    assert "cycle 3" in (output / "neb_profile.svg").read_text()
    candidate = json.loads((output / "provisional_candidates.jsonl").read_text())
    assert candidate["images"] == [3]
    assert candidate["source_sha256"] == report.source_sha256
    assert candidate["coordinate_association"] == "unverified; no release launched"
    archive = output / "reports" / "report_000001_cycle_000003"
    assert (archive / "neb_profile.svg").read_bytes() == (
        output / "neb_profile.svg"
    ).read_bytes()
    watcher.update_progress_once(source, output, 5)
    assert len((output / "provisional_candidates.jsonl").read_text().splitlines()) == 1
    source.write_bytes(source.read_bytes() + native_report([-10, -8, -9, -7, -11]))
    watcher.update_progress_once(source, output, 5)
    assert len((output / "provisional_candidates.jsonl").read_text().splitlines()) == 2


def test_once_and_flat_profile(tmp_path: Path) -> None:
    source = tmp_path / "neb.mdout"
    source.write_bytes(native_report([-10, -10, -10]))
    output = tmp_path / "watch"
    assert (
        watcher.main(
            [
                "--source",
                str(source),
                "--destination",
                str(output),
                "--expected-images",
                "3",
                "--once",
            ]
        )
        == 0
    )
    assert not (output / "provisional_candidates.jsonl").exists()
    assert "pending convergence marker" in (output / "neb_profile.svg").read_text()
