"""Tests for the active-site Amber measurement CLI."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import measure_nyl12_active_site as cli


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    topology = tmp_path / "system.parm7"
    coordinates = tmp_path / "system.rst7"
    topology.write_text("topology", encoding="utf-8")
    coordinates.write_text("coordinates", encoding="utf-8")
    return topology, coordinates


def _mock_measurement(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli, "AmberParm", lambda *args, **kwargs: object())
    monkeypatch.setattr(
        cli, "coordinates_from_amber", lambda model, target: {"xyz": (1, 2, 3)}
    )
    monkeypatch.setattr(
        cli,
        "measure_active_site",
        lambda coordinates, substrate, state: {
            "state": state,
            "substrate": substrate,
            "metrics": {"d1": 1.0},
            "geometry_source_sha256": "hash:reference_math.py",
        },
    )
    monkeypatch.setattr(cli, "_sha256", lambda path: f"hash:{path.name}")


def test_measure_files_writes_report_with_input_and_code_hashes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _mock_measurement(monkeypatch)
    topology, coordinates = _inputs(tmp_path)
    output = tmp_path / "reports" / "pa6.json"
    report = cli.measure_files(topology, coordinates, "pa6_l2", "ti1", output)
    saved = json.loads(output.read_text(encoding="utf-8"))
    assert saved == report
    assert report["input_hashes"] == {
        "topology_sha256": "hash:system.parm7",
        "coordinates_sha256": "hash:system.rst7",
    }
    assert report["input_paths"] == {
        "topology": str(topology.resolve()),
        "coordinates": str(coordinates.resolve()),
    }
    assert report["source_hashes"] == {
        "measurement_script_sha256": "hash:measure_nyl12_active_site.py",
        "active_site_geometry_sha256": "hash:active_site_geometry.py",
        "reference_math_sha256": "hash:reference_math.py",
    }
    assert report["target"] == "pa6_l2"


def test_measure_files_refuses_overwrite_and_invalid_selection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _mock_measurement(monkeypatch)
    topology, coordinates = _inputs(tmp_path)
    output = tmp_path / "report.json"
    cli.measure_files(topology, coordinates, "pa6_l2", "es", output)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        cli.measure_files(topology, coordinates, "pa6_l2", "es", output)
    with pytest.raises(ValueError, match="Unsupported target"):
        cli.measure_files(topology, coordinates, "pa6", "es", tmp_path / "other.json")
    with pytest.raises(ValueError, match="state"):
        cli.measure_files(
            topology, coordinates, "pa6_l2", "unknown", tmp_path / "other.json"
        )


def test_measure_files_requires_both_input_files(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        cli.measure_files(
            tmp_path / "missing.parm7",
            tmp_path / "missing.rst7",
            "pa6_l2",
            "ti1",
            tmp_path / "report.json",
        )


def test_parser_requires_explicit_state_and_output() -> None:
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["--topology", "a", "--coordinates", "b", "--target", "pa6_l2"]
        )
    parsed = parser.parse_args(
        [
            "--topology",
            "a",
            "--coordinates",
            "b",
            "--target",
            "poly_xama_maxa",
            "--state",
            "es",
            "--output",
            "o",
        ]
    )
    assert parsed.state == "es"


def test_main_uses_parsed_arguments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    called: dict[str, object] = {}
    monkeypatch.setattr(
        cli, "measure_files", lambda *args: called.setdefault("args", args)
    )
    topology = tmp_path / "in.parm7"
    coordinates = tmp_path / "in.rst7"
    output = tmp_path / "out.json"
    assert (
        cli.main(
            [
                "--topology",
                str(topology),
                "--coordinates",
                str(coordinates),
                "--target",
                "poly_xama_maxa",
                "--state",
                "ti1",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    assert called["args"] == (topology, coordinates, "poly_xama_maxa", "ti1", output)
