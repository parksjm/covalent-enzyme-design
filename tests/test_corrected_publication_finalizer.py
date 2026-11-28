"""Qualify exact curated transfer and publication-only pull exposure."""

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

SOURCE = (
    Path(__file__).resolve().parents[1]
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/boltz_masked_loop95_designs_20261127/template_alignment_repair01/finalize_publication.py"
)
spec = importlib.util.spec_from_file_location("repair_finalizer", SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fixture(tmp_path: Path) -> tuple[Path, Path]:
    output = tmp_path / "inspection_export_20261127"
    output.mkdir()
    (output / "query.yaml").write_text("sequences: []\n")
    digest = hashlib.sha256((output / "query.yaml").read_bytes()).hexdigest()
    (output / "MANIFEST.sha256").write_text(f"{digest}  query.yaml\n")
    (tmp_path / "publication_acceptance.json").write_text(
        json.dumps(
            {
                "status": "READY_TO_PULL",
                "matched_cases": 6,
                "target_phi_sign": 1,
                "manifest_sha256": hashlib.sha256(
                    (output / "MANIFEST.sha256").read_bytes()
                ).hexdigest(),
            }
        )
    )
    (tmp_path / "pull_after_acceptance.block").write_text(
        "# unique corrected block\necho accepted\n"
    )
    pull = tmp_path / "pull.sh"
    pull.write_text("#!/bin/bash\necho prior\n")
    return tmp_path, pull


def test_filtered_yaml_transfer_and_append(tmp_path: Path) -> None:
    base, pull = fixture(tmp_path)
    m.verify_and_wire(base, pull)
    assert pull.read_text().startswith("#!/bin/bash\necho prior\n")
    assert (
        json.loads((base / "filtered_pull_acceptance.json").read_text())["status"]
        == "accepted"
    )
    with pytest.raises(FileExistsError):
        m.verify_and_wire(base, pull)


def test_reject_incomplete(tmp_path: Path) -> None:
    base, pull = fixture(tmp_path)
    (base / "publication_acceptance.json").write_text("{}")
    with pytest.raises(ValueError):
        m.verify_and_wire(base, pull)
    assert "corrected block" not in pull.read_text()


def test_fail_filtered_unpublished_file(tmp_path: Path) -> None:
    base, pull = fixture(tmp_path)
    output = base / "inspection_export_20261127"
    (output / "query.yaml").write_text("corrupted")
    with pytest.raises(__import__("subprocess").CalledProcessError):
        m.verify_and_wire(base, pull)
    assert "corrected block" not in pull.read_text()


def test_main_binding(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    base, pull = fixture(tmp_path)
    source = base / "publish_corrected_comparison.py"
    source.write_text("pass\n")
    binding = base / "followthrough_source_binding.json"
    binding.write_text(
        json.dumps(
            {"publisher_sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
        )
    )
    monkeypatch.setattr(m, "BASE", base)
    monkeypatch.setattr(m, "ROOT", base)
    calls = []
    monkeypatch.setattr(m.subprocess, "run", lambda *a, **k: calls.append(a))
    monkeypatch.setattr(m, "verify_and_wire", lambda *a: calls.append(a))
    m.main()
    assert len(calls) == 2
    source.write_text("changed")
    with pytest.raises(ValueError):
        m.main()
