"""Tests for packaging explicitly accepted relay-repair outputs."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pytest

from scripts import package_nyl12_relay_inspection as package_cli


def _result(path: Path, target: str, *, accepted: bool = True) -> None:
    """Write the minimal acceptance and target metadata consumed by packaging."""
    path.mkdir(parents=True)
    (path / "acceptance.json").write_text(
        json.dumps({"accepted": accepted}), encoding="utf-8"
    )
    (path / f"{target}_ti1_provenance.json").write_text(
        json.dumps({"target": target}), encoding="utf-8"
    )


def _patch_package_inputs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> tuple[dict[str, Any], dict[str, Path]]:
    """Patch repository selection and package build, capturing selected sources."""
    selected = {"existing/retained.txt": tmp_path / "retained.txt"}
    selected["existing/retained.txt"].write_text("retained", encoding="utf-8")
    captured: dict[str, Any] = {}

    def selected_files(*args: Any) -> dict[str, Path]:
        captured["selected_args"] = args
        return dict(selected)

    def build_package(files: dict[str, Path], output: Path) -> dict[str, Any]:
        captured["files"] = files
        captured["output"] = output
        return {"files": sorted(files), "payload_bytes": 17}

    monkeypatch.setattr(package_cli, "selected_files", selected_files)
    monkeypatch.setattr(package_cli, "build_package", build_package)
    return captured, selected


def test_package_repairs_includes_only_explicitly_accepted_results(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Accepted target reports and reproducibility sources reach package build."""
    captured, selected = _patch_package_inputs(monkeypatch, tmp_path)
    repo = tmp_path / "repo"
    pa6 = tmp_path / "pa6"
    maxa = tmp_path / "maxa"
    journal = tmp_path / "actions.jsonl"
    journal.write_text("journal", encoding="utf-8")
    _result(pa6, "pa6_l4")
    _result(maxa, "maxa_maxa")
    output = tmp_path / "bundle"

    result = package_cli.package_repairs(repo, pa6, maxa, journal, output)

    files = captured["files"]
    assert result == {"files": sorted(files), "payload_bytes": 17}
    assert files["existing/retained.txt"] == selected["existing/retained.txt"]
    assert files["provenance/agent_actions_at_export.jsonl"] == journal
    for target, directory in (("pa6_l4", pa6), ("maxa_maxa", maxa)):
        assert files[f"protein/{target}_ti1_starting.pdb"] == directory / f"{target}_ti1_starting.pdb"
        assert files[f"reports/{target}_relay_repair/acceptance.json"] == directory / "acceptance.json"
        assert files[f"reports/{target}_original_active_site.json"] == directory / "original_active_site.json"
        assert files[f"provenance/{target}_repair_job/status.json"] == directory.parent / "job/status.json"
    assert captured["output"] == output


@pytest.mark.parametrize(
    ("bad_target", "accepted", "metadata_target", "message"),
    [
        ("pa6_l4", False, "pa6_l4", "unaccepted repair for pa6_l4"),
        ("maxa_maxa", True, "wrong_target", "Repair target mismatch"),
    ],
)
def test_package_repairs_rejects_unaccepted_or_mismatched_result(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    bad_target: str,
    accepted: bool,
    metadata_target: str,
    message: str,
) -> None:
    """Failed acceptance and provenance-target mismatch stop before packaging."""
    captured, _ = _patch_package_inputs(monkeypatch, tmp_path)
    pa6 = tmp_path / "pa6"
    maxa = tmp_path / "maxa"
    _result(pa6, "pa6_l4", accepted=accepted if bad_target == "pa6_l4" else True)
    _result(
        maxa,
        metadata_target if bad_target == "maxa_maxa" else "maxa_maxa",
        accepted=accepted if bad_target == "maxa_maxa" else True,
    )
    if bad_target == "maxa_maxa" and metadata_target != "maxa_maxa":
        (maxa / "wrong_target_ti1_provenance.json").unlink()
        (maxa / "maxa_maxa_ti1_provenance.json").write_text(
            json.dumps({"target": metadata_target}), encoding="utf-8"
        )
    journal = tmp_path / "actions.jsonl"
    journal.write_text("journal", encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        package_cli.package_repairs(
            tmp_path / "repo", pa6, maxa, journal, tmp_path / "bundle"
        )

    assert "files" not in captured


def test_main_forwards_arguments_to_package_repairs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The CLI forwards all explicitly selected paths to its package builder."""
    arguments = ["repo", "pa6", "maxa", "journal", "bundle"]
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "package_nyl12_relay_inspection.py",
            "--repo", arguments[0],
            "--pa6-result", arguments[1],
            "--maxa-result", arguments[2],
            "--journal", arguments[3],
            "--output", arguments[4],
        ],
    )
    captured: dict[str, Any] = {}

    def package_repairs(*paths: Path) -> dict[str, Any]:
        captured["paths"] = paths
        return {"files": ["one"], "payload_bytes": 17}

    monkeypatch.setattr(package_cli, "package_repairs", package_repairs)
    package_cli.main()

    assert captured["paths"] == tuple(Path(value) for value in arguments)
