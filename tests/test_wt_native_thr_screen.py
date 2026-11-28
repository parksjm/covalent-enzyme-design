"""Tests for compact full-eligible-frame native-Thr screening."""

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from nyl12_covalent import wt_native_thr_screen as screen_module


def test_eligible_frames_retains_native_failures_and_rejects_bad_time(
    tmp_path: Path,
) -> None:
    path = tmp_path / "screen.csv"
    path.write_text(
        "frame,time_ns,tier,native_torsions,psi_delta,chi1_delta\n"
        "1.0,0.05,A,True,-10,5\n2.0,0.10,B,False,-50,5\n"
        "3.0,0.15,none,True,-10,5\n4.0,0.20,B,True,20,-10\n",
        encoding="utf-8",
    )
    selected = screen_module.eligible_frames(path, "r1")
    assert [(row["frame"], row["tier"]) for row in selected] == [
        (1, "A"),
        (2, "B"),
        (4, "B"),
    ]
    assert selected[1]["native_torsions"] is False
    path.write_text(
        "frame,time_ns,tier,native_torsions,psi_delta,chi1_delta\n1,nan,A,True,-10,5\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Invalid time"):
        screen_module.eligible_frames(path, "r1")


def test_cpptraj_text_extracts_only_listed_frames(tmp_path: Path) -> None:
    topology = tmp_path / "system.prmtop"
    trajectory = tmp_path / "production.nc"
    output = tmp_path / "eligible.nc"
    frames = [
        {"replica": "r1", "frame": 2},
        {"replica": "r1", "frame": 19},
    ]
    text = screen_module.cpptraj_text(
        topology, frames, {"r1": trajectory}, ":1313-1316", output
    )
    assert f"trajin {trajectory.resolve()} 2 2 1" in text
    assert f"trajin {trajectory.resolve()} 19 19 1" in text
    assert "trajout" in text and "autoimage anchor :1313-1316" in text
    with pytest.raises(ValueError, match="anchor"):
        screen_module.cpptraj_text(topology, frames, {"r1": trajectory}, "", output)


def test_candidate_rank_ignores_native_similarity() -> None:
    def candidate(
        og1: float, rmsd: float, n_shift: float, frame: int
    ) -> dict[str, Any]:
        return {
            "tier": "B",
            "time_ns": float(frame),
            "frame": frame,
            "native_alignment": {
                "thr233_held_out": {
                    "heavy_atom_displacements_angstrom": {"OG1": og1, "N": n_shift},
                    "rmsd_angstrom": rmsd,
                }
            },
        }

    ordered = sorted(
        [
            candidate(0.1, 0.1, 0.1, 3),
            candidate(0.5, 1.0, 0.7, 2),
            candidate(4.0, 4.0, 4.0, 1),
        ],
        key=screen_module._candidate_rank,
    )
    assert [item["frame"] for item in ordered] == [1, 2, 3]
    assert ordered[0]["native_alignment"]["thr233_held_out"]["rmsd_angstrom"] == 4.0


def test_alignment_helpers_load_mapping_once_and_measure_frame(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference_path = tmp_path / "reference.cif"
    topology_path = tmp_path / "topology.prmtop"
    reference_path.write_text("reference", encoding="utf-8")
    topology_path.write_text("topology", encoding="utf-8")
    mapped_reference = {
        (auth, atom): ("C", np.array([float(auth), 0.0, 0.0]))
        for _, auth, _ in screen_module.SITE_RESIDUES
        for atom in screen_module.BACKBONE
    }
    monkeypatch.setattr(screen_module, "_cif_atoms", lambda *args: mapped_reference)
    monkeypatch.setattr(screen_module.parmed, "load_file", lambda _: SimpleNamespace())
    reference, structure, fit_keys, reference_fit = screen_module._alignment_inputs(
        reference_path, topology_path
    )
    assert reference is mapped_reference
    assert structure is not None
    assert len(fit_keys) == 16
    assert reference_fit.shape == (16, 3)

    structure = SimpleNamespace(atoms=[object()], coordinates=None)
    mobile = {key: value for key, value in mapped_reference.items()}
    monkeypatch.setattr(screen_module, "_topology_atoms", lambda *_: mobile)
    monkeypatch.setattr(
        screen_module,
        "fit_coordinates",
        lambda *_: (np.eye(3), np.zeros(3), np.zeros(3), 0.25),
    )
    expected = {"thr233_held_out": {"rmsd_angstrom": 0.5}}
    monkeypatch.setattr(screen_module, "state_metrics", lambda *_: expected)
    result = screen_module._measure_frame(
        np.array([[1.0, 2.0, 3.0]]),
        reference,
        structure,
        fit_keys,
        reference_fit,
    )
    assert result == {"fit_rmsd_angstrom": 0.25, **expected}
    with pytest.raises(ValueError, match="do not match"):
        screen_module._measure_frame(
            np.zeros((2, 3)), reference, structure, fit_keys, reference_fit
        )


def test_screen_combines_candidates_and_writes_ranked_extraction_list(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source_paths = {}
    for name in ("topology", "reference", "cpptraj", "trajectory", "screen"):
        path = tmp_path / name
        path.write_text("source", encoding="utf-8")
        source_paths[name] = path
    source_paths["screen"].write_text(
        "frame,time_ns,tier,native_torsions,psi_delta,chi1_delta\n"
        "1,0.05,A,True,-10,5\n2,0.10,B,False,-50,5\n3,0.15,B,True,20,-10\n",
        encoding="utf-8",
    )
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "topology": str(source_paths["topology"]),
                "reference_cif": str(source_paths["reference"]),
                "image_anchor": ":1313-1316",
                "replicas": [
                    {
                        "name": "r1",
                        "trajectory": str(source_paths["trajectory"]),
                        "screen_csv": str(source_paths["screen"]),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    def fake_run(command: list[str], **kwargs: Any) -> None:
        input_text = Path(command[-1]).read_text(encoding="utf-8")
        assert "trajin" in input_text
        Path(input_text.split("trajout ", 1)[1].split()[0]).write_bytes(b"compact")

    class FakeTrajectory:
        coordinates = np.array(
            [[[0.1, 0.0, 0.0]], [[0.2, 0.0, 0.0]], [[0.3, 0.0, 0.0]]]
        )

        def close(self) -> None:
            return None

    monkeypatch.setattr(screen_module.subprocess, "run", fake_run)
    monkeypatch.setattr(
        screen_module.NetCDFTraj, "open_old", lambda _: FakeTrajectory()
    )
    monkeypatch.setattr(
        screen_module,
        "_alignment_inputs",
        lambda *args: ({}, SimpleNamespace(), [], np.empty((0, 3))),
    )

    def measure(coordinates: np.ndarray, *_: Any) -> dict[str, Any]:
        displacement = float(coordinates[0, 0])
        n_displacement = displacement + 0.1 if displacement < 0.15 else 1.3
        og1_displacement = displacement if displacement < 0.15 else 0.9
        return {
            "fit_rmsd_angstrom": 0.1,
            "thr233_held_out": {
                "heavy_atom_displacements_angstrom": {
                    "OG1": og1_displacement,
                    "N": n_displacement,
                },
                "rmsd_angstrom": displacement + 0.2,
            },
        }

    monkeypatch.setattr(screen_module, "_measure_frame", measure)
    output = tmp_path / "out"
    report = screen_module.screen(config, output, source_paths["cpptraj"], 3)
    assert report["eligible_frame_count"] == 3
    assert report["accepted_count"] == 3
    assert report["native_similarity_count"] == 1
    assert report["selection_policy"]["native_thr_role"] == "diagnostic_only"
    assert report["best_candidates_for_extraction"][0]["frame"] == 1
    assert report["ranked_candidates"][1]["accepted_native_gate"] is False
    assert report["ranked_candidates"][1]["selection_eligible"] is True
    assert [r["frame"] for r in report["best_candidates_for_extraction"]] == [1, 2, 3]
    assert (output / "eligible_candidates.nc").read_bytes() == b"compact"
    assert (output / "ranked_candidates.csv").is_file()
    assert json.loads((output / "native_thr_screen.json").read_text())["top_k"] == 3
