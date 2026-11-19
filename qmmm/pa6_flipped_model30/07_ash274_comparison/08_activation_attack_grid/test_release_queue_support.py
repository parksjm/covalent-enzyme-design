"""Focused checks for retained-chemistry candidate selection."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import release_queue_support as support


def geometry_values(**changes: float) -> dict[str, float]:
    """Return one valid chemistry record with optional coordinate changes."""
    values = {
        "q_N_A": -0.7,
        "attack_A": 2.0,
        "thr_n_hg1_A": 2.0,
        "thr_n_h1_A": 1.02,
        "thr_n_h2_A": 1.02,
        "ash274_od2_hd2_A": 1.0,
        "ash274_hd2_water_o_A": 5.0,
        "water_o_h1_A": 0.97,
        "leaving_n_hn_A": 1.02,
        "substrate_c_o_A": 1.26,
        "substrate_c_n_A": 1.36,
    }
    values.update(changes)
    return values


def write_candidate(
    root: Path, index: int, geometry: dict[str, float], energy: float
) -> dict[str, Any]:
    """Write a minimal hashed point/process/acceptance fixture."""
    point_id = f"p{index:02d}"
    restart = root / f"{point_id}.ncrst"
    restart.write_bytes(point_id.encode())
    restart_hash = hashlib.sha256(restart.read_bytes()).hexdigest()
    process_path = root / f"{point_id}_process.json"
    process_path.write_text(
        json.dumps(
            {
                "status": "finished",
                "exit_status": 0,
                "force_converged": True,
                "scc_converged": True,
            }
        )
    )
    acceptance_path = root / f"{point_id}_acceptance.json"
    acceptance_path.write_text(
        json.dumps(
            {
                "accepted": True,
                "force_converged": True,
                "scc_converged": True,
                "geometry": geometry,
                "energies": {"physical_energy_kcal_mol": energy},
                "restart_sha256": restart_hash,
            }
        )
    )
    return {
        "status": "completed",
        "restart": restart.name,
        "restart_sha256": restart_hash,
        "process_record": process_path.name,
        "acceptance_record": acceptance_path.name,
        "chemistry_class": "retained_candidate",
    }


def test_inventory_requires_thr_and_ash_proton_retention() -> None:
    """Reject lost Thr N or ASH274 protons from release selection."""
    assert support.inventory_retained(geometry_values())
    assert not support.inventory_retained(geometry_values(thr_n_h2_A=1.40))
    assert not support.inventory_retained(geometry_values(ash274_od2_hd2_A=1.40))


def test_candidate_classes_use_independent_geometry_gates() -> None:
    """Recognize only the requested activated and tetrahedral-like geometries."""
    activated = geometry_values(q_N_A=0.51, thr_n_hg1_A=1.20)
    ti_like = geometry_values(
        thr_n_hg1_A=1.20,
        attack_A=1.70,
        substrate_c_o_A=1.31,
        substrate_c_n_A=1.54,
    )
    assert support.candidate_classes(activated) == ["activated"]
    assert support.candidate_classes(ti_like) == ["ti_like"]
    assert support.candidate_classes(geometry_values(q_N_A=0.5)) == []


def test_selects_lowest_energy_candidate_per_unique_class(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Choose each class minimum and omit an absent class explicitly."""
    points = {
        "points": {
            "activated_high": write_candidate(
                tmp_path, 1, geometry_values(q_N_A=0.7, thr_n_hg1_A=1.0), -10.0
            ),
            "activated_low": write_candidate(
                tmp_path, 2, geometry_values(q_N_A=0.8, thr_n_hg1_A=1.0), -20.0
            ),
        }
    }
    for index in range(3, 26):
        points["points"][f"unused_{index}"] = {"status": "failed"}
    point_path = tmp_path / "grid_points.json"
    point_path.write_text(json.dumps(points))
    monkeypatch.setattr(support, "GRID", tmp_path)
    monkeypatch.setattr(support, "POINTS", point_path)

    selected, decisions = support.select_candidates()

    assert [item["key"] for item in selected] == ["activated_low"]
    assert decisions["classes"]["activated"]["status"] == "selected"
    assert decisions["classes"]["ti_like"]["status"] == "omitted"
