"""Tests for reaction path measurements and inspection exports."""

from __future__ import annotations

import json
import struct
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import path_inspection as inspection
from nyl12_covalent.path_distances import distance_jumps, measure_distances


def _write_dcd(path: Path, frames: int = 1, atoms: int = 3) -> None:
    header = b"CORD" + struct.pack("<i", frames) + bytes(72) + struct.pack("<i", 0)
    with path.open("wb") as stream:
        stream.write(struct.pack("<i", 84))
        stream.write(header)
        stream.write(struct.pack("<i", 84))
        stream.write(struct.pack("<i", 4))
        stream.write(struct.pack("<i", 4))
        stream.write(struct.pack("<i", 4))
        stream.write(struct.pack("<i", 4))
        stream.write(struct.pack("<i", atoms))
        stream.write(struct.pack("<i", 4))


def test_measure_distances_uses_one_based_indices() -> None:
    coordinates = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]])
    assert measure_distances(coordinates, {"bond": {"atoms": [1, 2]}}) == {
        "bond": 5.0
    }


@pytest.mark.parametrize(
    ("coordinates", "pairs"),
    [
        (np.zeros((2, 2)), {"bond": {"atoms": [1, 2]}}),
        (np.array([[0.0, np.nan, 0.0]]), {"bond": {"atoms": [1, 1]}}),
        (np.zeros((2, 3)), {}),
        (np.zeros((2, 3)), {"bond": {"atoms": [0, 2]}}),
        (np.zeros((2, 3)), {"bond": {"atoms": [1, 3]}}),
        (np.zeros((2, 3)), {"bond": {"atoms": [2, 2]}}),
    ],
)
def test_measure_distances_rejects_invalid_inputs(
    coordinates: np.ndarray, pairs: dict[str, dict[str, object]]
) -> None:
    with pytest.raises(ValueError):
        measure_distances(coordinates, pairs)  # type: ignore[arg-type]


def test_distance_jumps_keeps_stage_boundary_and_signed_change() -> None:
    rows = [
        {"window": "a", "stage": "direct", "distances_A": {"bond": 1.0}},
        {"window": "b", "stage": "relay", "distances_A": {"bond": 1.4}},
    ]
    assert distance_jumps(rows, 0.3) == [
        {
            "from": "a",
            "to": "b",
            "distance": "bond",
            "delta_A": pytest.approx(0.4),
            "stage_boundary": True,
        }
    ]


@pytest.mark.parametrize("threshold", [0.0, -1.0, float("nan"), float("inf")])
def test_distance_jumps_rejects_invalid_threshold(threshold: float) -> None:
    with pytest.raises(ValueError, match="threshold"):
        distance_jumps([], threshold)


def test_distance_jumps_rejects_changed_pairs_and_nonfinite_delta() -> None:
    first = {"window": "a", "stage": "direct", "distances_A": {"a": 1.0}}
    with pytest.raises(ValueError, match="definitions changed"):
        distance_jumps([first, {**first, "window": "b", "distances_A": {"b": 2.0}}], 1)
    with pytest.raises(ValueError, match="Nonfinite"):
        distance_jumps([first, {**first, "window": "b", "distances_A": {"a": np.nan}}], 1)


def test_dcd_counts_accepts_standard_coordinate_header(tmp_path: Path) -> None:
    dcd = tmp_path / "path.dcd"
    _write_dcd(dcd, frames=4, atoms=19)
    assert inspection.dcd_counts(dcd) == (4, 19)


@pytest.mark.parametrize("payload", [b"bad", struct.pack("<i", 84) + b"CORD"])
def test_dcd_counts_rejects_truncated_or_invalid_header(
    tmp_path: Path, payload: bytes
) -> None:
    dcd = tmp_path / "bad.dcd"
    dcd.write_bytes(payload)
    with pytest.raises(ValueError):
        inspection.dcd_counts(dcd)


def _write_endpoint(run: Path, label: str, stage: str, converged: bool) -> None:
    folder = run / label
    folder.mkdir()
    (folder / "optimized.ncrst").write_bytes(b"coordinates")
    (folder / "acceptance.json").write_text(
        json.dumps({"native_converged": converged})
    )
    (folder / "process.json").write_text(
        json.dumps({"exit_status": 0, "stage": stage})
    )


def test_collect_endpoints_includes_accepted_and_explicitly_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = tmp_path / "run"
    run.mkdir()
    topology = tmp_path / "system.parm7"
    topology.write_bytes(b"topology")
    (run / "status.json").write_text(
        json.dumps({"state": "stopped", "completed_windows": ["w1"]})
    )
    (run / "inputs.json").write_text(
        json.dumps({"topology": str(topology), "box": [10, 10, 10, 90, 90, 90]})
    )
    _write_endpoint(run, "w1", "direct", True)
    _write_endpoint(run, "w2", "relay", False)

    def fake_parm(*args: str) -> SimpleNamespace:
        return SimpleNamespace(
            atoms=[
                SimpleNamespace(name="C1", residue=SimpleNamespace(idx=0)),
                SimpleNamespace(name="O1", residue=SimpleNamespace(idx=0)),
            ],
            coordinates=np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]),
        )

    monkeypatch.setattr(inspection, "AmberParm", fake_parm)
    monkeypatch.setattr(
        inspection,
        "Rst7",
        SimpleNamespace(
            open=lambda restart: SimpleNamespace(
                coordinates=np.array([[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]])
            )
        ),
        raising=False,
    )
    config = {
        "include_rejected_endpoint": True,
        "atom_identities": {
            "1": {"name": "C1", "residue": 1},
            "2": {"name": "O1", "residue": 1},
        },
        "pairs": {"bond": {"atoms": [1, 2]}},
        "jump_threshold_A": 0.2,
    }
    result = inspection.collect_endpoints(run, config)
    assert [(row["window"], row["accepted"]) for row in result["rows"]] == [
        ("w1", True),
        ("w2", False),
    ]
    assert result["rows"][0]["distances_A"] == {"bond": pytest.approx(1.2)}


def test_collect_endpoints_rejects_failed_engine_and_identity_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = tmp_path / "run"
    run.mkdir()
    topology = tmp_path / "system.parm7"
    topology.write_bytes(b"topology")
    (run / "status.json").write_text(json.dumps({"state": "done", "completed_windows": ["w1"]}))
    (run / "inputs.json").write_text(json.dumps({"topology": str(topology), "box": [1] * 6}))
    _write_endpoint(run, "w1", "direct", True)
    config = {
        "include_rejected_endpoint": False,
        "atom_identities": {"1": {"name": "wrong", "residue": 1}},
        "pairs": {"bond": {"atoms": [1, 2]}},
        "jump_threshold_A": 0.2,
    }
    monkeypatch.setattr(
        inspection,
        "AmberParm",
        lambda *_args: SimpleNamespace(
            atoms=[SimpleNamespace(name="C1", residue=SimpleNamespace(idx=0))],
            coordinates=np.zeros((1, 3)),
        ),
    )
    monkeypatch.setattr(
        inspection,
        "Rst7",
        SimpleNamespace(open=lambda _restart: SimpleNamespace(coordinates=np.zeros((1, 1, 3)))),
        raising=False,
    )
    with pytest.raises(ValueError, match="identity mismatch"):
        inspection.collect_endpoints(run, config)
    config["atom_identities"] = {"1": {"name": "C1", "residue": 1}}
    (run / "w1" / "process.json").write_text(json.dumps({"exit_status": 1}))
    with pytest.raises(ValueError, match="failed engine"):
        inspection.collect_endpoints(run, config)


def test_publish_inspection_validates_and_publishes_mocked_exports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run = tmp_path / "run"
    run.mkdir()
    (run / "inputs.json").write_text(json.dumps({"topology": "system.parm7"}))
    config_path = tmp_path / "config.json"
    plot_script = tmp_path / "plot.py"
    plot_script.write_text("# plot worker")
    assets = tmp_path / "assets"
    assets.mkdir()
    for name in ("stripped.parm7.gz", "qm_indices.tcl", "nyl12_es_vmd_common.tcl"):
        (assets / name).write_text(name)
    output = tmp_path / "inspection"
    config_path.write_text(
        json.dumps(
            {
                "output": str(output),
                "pairs": {"bond": {"atoms": [1, 2]}},
                "viewer_assets": str(assets),
                "export_environment": {"OMP_NUM_THREADS": "1"},
                "cpptraj": "cpptraj",
                "plot_python": "python",
                "strip_mask": ":WAT",
                "expected_stripped_atoms": 3,
            }
        )
    )
    rows = [
        {
            "window": "w1",
            "stage": "direct",
            "accepted": True,
            "restart": "w1/optimized.ncrst",
            "distances_A": {"bond": 1.3},
        }
    ]
    monkeypatch.setattr(
        inspection,
        "collect_endpoints",
        lambda *_: {
            "rows": rows,
            "source_sha256": {},
            "run_state": {"state": "done"},
            "pairs": {"bond": {"atoms": [1, 2]}},
            "box": [10, 10, 10, 90, 90, 90],
            "jump_threshold_A": 0.3,
            "jumps": [],
        },
    )

    def fake_run(command: list[str], **kwargs: object) -> None:
        cwd = Path(kwargs["cwd"]) if "cwd" in kwargs else None
        if cwd is not None:
            _write_dcd(cwd / "reaction_path.dcd", frames=1, atoms=3)
        else:
            args = command
            Path(args[-1]).with_name("reaction_distances.png").write_bytes(b"PNG")

    monkeypatch.setattr(inspection.subprocess, "run", fake_run)
    generation = inspection.publish_inspection(run, config_path, plot_script)
    assert generation.name == "update_0001"
    assert (output / "latest").resolve() == generation
    assert inspection.dcd_counts(generation / "reaction_path.dcd") == (1, 3)
    assert (generation / "reaction_distances.csv").is_file()
    assert (generation / "reaction_distances.png").read_bytes() == b"PNG"
    assert (generation / "MANIFEST.sha256").is_file()
