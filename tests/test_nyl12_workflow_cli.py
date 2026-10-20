from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts import audit_nyl12_ti1_qm_region as audit_cli
from scripts import run_nyl12_job as job_cli
from scripts import search_nyl12_relay_water as water_cli


def test_qm_audit_writes_json_and_passes_provenance_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    report = {"qmcharge": 0, "topology_sha256": "abc", "sources": ["source"]}
    received: list[object] = []
    monkeypatch.setattr(
        audit_cli, "audit_qm_region", lambda *args: received.extend(args) or report
    )
    output = tmp_path / "audit.json"
    monkeypatch.setattr(
        "sys.argv",
        [
            "audit_nyl12_ti1_qm_region.py",
            "input.prmtop",
            "--source",
            "input.inpcrd",
            "--source",
            "input.json",
            "--output",
            str(output),
        ],
    )

    audit_cli.main()

    assert received == [
        Path("input.prmtop"),
        [Path("input.inpcrd"), Path("input.json")],
    ]
    assert json.loads(output.read_text()) == report


def test_qm_audit_prints_json_when_output_is_omitted(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(audit_cli, "audit_qm_region", lambda *_: {"qmcharge": 0})
    monkeypatch.setattr("sys.argv", ["audit_nyl12_ti1_qm_region.py", "input.prmtop"])

    audit_cli.main()

    assert json.loads(capsys.readouterr().out) == {"qmcharge": 0}


def test_job_status_prints_existing_status_json(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    output = tmp_path / "job"
    output.mkdir()
    status = {"session_name": "audit", "exit_status": 0}
    (output / "status.json").write_text(json.dumps(status))

    assert job_cli.main(["status", str(output)]) == 0
    assert json.loads(capsys.readouterr().out) == status


def test_job_launch_parses_environment_inputs_and_separator(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_launch(*args: object, **kwargs: object) -> Path:
        captured["args"] = args
        captured["kwargs"] = kwargs
        return tmp_path / "run.sh"

    monkeypatch.setattr(job_cli, "launch_job", fake_launch)
    result = job_cli.main(
        [
            "launch",
            "--session",
            "unit_job",
            "--output-dir",
            str(tmp_path / "out"),
            "--cwd",
            str(tmp_path),
            "--env",
            "OMP_NUM_THREADS=4",
            "--env",
            "LABEL=a=b",
            "--input",
            "input.pdb",
            "--",
            "python",
            "worker.py",
        ]
    )

    assert result == 0
    assert captured["args"] == (
        "unit_job",
        tmp_path / "out",
        tmp_path,
        ["python", "worker.py"],
    )
    assert captured["kwargs"] == {
        "env_overrides": {"OMP_NUM_THREADS": "4", "LABEL": "a=b"},
        "input_paths": [Path("input.pdb")],
    }


def test_job_launch_rejects_environment_without_equals(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="KEY=VALUE"):
        job_cli.main(
            [
                "launch",
                "--session",
                "unit_job",
                "--output-dir",
                str(tmp_path / "out"),
                "--cwd",
                str(tmp_path),
                "--env",
                "MISSING_EQUALS",
                "--",
                "python",
            ]
        )


class FakeTopology:
    def __init__(self, *_: object, **__: object) -> None:
        self.coordinates = np.zeros((4, 3), dtype=float)
        self.saved_coordinates: list[np.ndarray] = []

    def save(self, path: str, format: str) -> None:
        assert format == "rst7"
        self.saved_coordinates.append(self.coordinates.copy())
        Path(path).write_text("synthetic restart\n")


def _water_report(accepted: bool) -> dict[str, object]:
    return {
        "water_atom_indices_zero_based": [1, 2],
        "accepted_oriented_candidate_count": int(accepted),
        "candidates": [
            {
                "accepted": accepted,
                "water_xyz_A": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                "maximum_contact_shortfall_A": 0.0 if accepted else 0.5,
            }
        ],
    }


def _run_water_search(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    report: dict[str, object],
) -> tuple[FakeTopology, Path]:
    topology_path = tmp_path / "input.prmtop"
    coordinates_path = tmp_path / "input.inpcrd"
    topology_path.write_text("topology\n")
    coordinates_path.write_text("coordinates\n")
    output = tmp_path / "search"
    topology = FakeTopology()
    monkeypatch.setattr(water_cli, "AmberParm", lambda *args, **kwargs: topology)
    monkeypatch.setattr(water_cli, "search_relay_water", lambda *_: report)
    monkeypatch.setattr(
        "sys.argv",
        [
            "search_nyl12_relay_water.py",
            "--topology",
            str(topology_path),
            "--coordinates",
            str(coordinates_path),
            "--oxyanion",
            "O49",
            "--output",
            str(output),
        ],
    )
    return topology, output


def test_water_search_saves_accepted_coordinates_and_hash_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    topology, output = _run_water_search(tmp_path, monkeypatch, _water_report(True))

    water_cli.main()

    saved = output / "candidate_01.rst7"
    assert saved.is_file()
    assert topology.saved_coordinates[0].tolist() == [
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [0.0, 0.0, 0.0],
    ]
    report = json.loads((output / "water_search.json").read_text())
    assert report["input_sha256"]
    assert str(Path(__file__).resolve()) not in report["input_sha256"]


def test_water_search_preserves_failed_report_then_raises(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, output = _run_water_search(tmp_path, monkeypatch, _water_report(False))

    with pytest.raises(RuntimeError, match="No relay-water seed"):
        water_cli.main()

    report = json.loads((output / "water_search.json").read_text())
    assert report["candidates"][0]["accepted"] is False
    assert report["input_sha256"]


def test_water_search_rejects_existing_output_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, output = _run_water_search(tmp_path, monkeypatch, _water_report(True))
    output.mkdir()

    with pytest.raises(FileExistsError):
        water_cli.main()
