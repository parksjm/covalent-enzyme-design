"""Test source ownership mapping and lossless restart extraction on tiny fixtures."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from scipy.io import netcdf_file
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/screen_nalpha_leaving_water.py"
SPEC = importlib.util.spec_from_file_location("nalpha_screen", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SCREEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCREEN)


def small_topology() -> Any:
    """Build an NH2 terminal N and one topology-owned three-atom water."""
    structure = parmed.Structure()
    for name, atomic_number, residue, number in [
        ("N", 7, "THR", 1),
        ("H1", 1, "THR", 1),
        ("H2", 1, "THR", 1),
        ("O", 8, "WAT", 2),
        ("H1", 1, "WAT", 2),
        ("H2", 1, "WAT", 2),
    ]:
        structure.add_atom(
            parmed.Atom(name=name, atomic_number=atomic_number), residue, number
        )
    for first, second in [(0, 1), (0, 2), (3, 4), (3, 5)]:
        structure.bonds.append(
            parmed.Bond(structure.atoms[first], structure.atoms[second])
        )
    return structure


def test_mapping_uses_bonded_hydrogen_owners_and_rejects_bad_state() -> None:
    """Residue-owned waters and N-H inventory must match the declared state."""
    structure = small_topology()
    config = {"atoms": {"nalpha": [1, "N"]}, "residue_names": {"1": "THR"}}
    mapped = SCREEN.mapping(structure, config)
    assert mapped["nalpha_h"] == [1, 2]
    assert mapped["waters"] == [(2, 3, 4, 5)]
    structure.atoms[2].atomic_number = 6
    with pytest.raises(ValueError, match="NH2"):
        SCREEN.mapping(structure, config)
    structure.atoms[2].atomic_number = 1
    structure.residues[0].name = "LYS"
    with pytest.raises(ValueError, match="Unexpected residue"):
        SCREEN.mapping(structure, config)
    with pytest.raises(ValueError, match="Ambiguous mapping"):
        SCREEN.mapping(structure, {"atoms": {"nalpha": [1, "NZ"]}})


def test_extract_preserves_coordinates_and_box_exactly(
    tmp_path: Path, monkeypatch: Any
) -> None:
    """NetCDF extraction must retain original floats, box, and ownership mapping."""
    structure = small_topology()
    topology = tmp_path / "topology.prmtop"
    topology.write_text("Synthetic fixture topology is monkeypatched in this test.\n")
    trajectory = tmp_path / "trajectory.nc"
    coordinates = np.arange(18, dtype=np.float32).reshape(1, 6, 3) / 7
    lengths = np.array([[20.0, 21.0, 22.0]])
    angles = np.array([[90.0, 100.0, 110.0]])
    with netcdf_file(str(trajectory), "w") as data:
        for dimension, count in [("frame", 1), ("atom", 6), ("xyz", 3)]:
            data.createDimension(dimension, count)
        data.createVariable("coordinates", "f", ("frame", "atom", "xyz"))[:] = (
            coordinates
        )
        data.createVariable("cell_lengths", "d", ("frame", "xyz"))[:] = lengths
        data.createVariable("cell_angles", "d", ("frame", "xyz"))[:] = angles
    monkeypatch.setattr(SCREEN.parmed, "load_file", lambda _: structure)
    config = {
        "cases": [
            {
                "name": "fixture",
                "topology": str(topology),
                "residue_names": {"1": "THR"},
            }
        ],
        "atoms": {"nalpha": [1, "N"]},
        "residue_names": {"1": "THR"},
    }
    row = {
        "case": "fixture",
        "frame": 1,
        "local_frame": 1,
        "trajectory": str(trajectory),
        "water_resid": 2,
    }
    record = SCREEN.extract(row, config, tmp_path, "selected")
    assert record["coordinates_unmodified"]
    assert record["atom_mapping"]["nalpha"]["index1"] == 1
    assert record["relay_water_atoms_index1"] == [4, 5, 6]
    restart = parmed.amber.NetCDFRestart.open_old(record["restart"])
    try:
        np.testing.assert_array_equal(
            np.asarray(restart.coordinates).reshape(6, 3), coordinates[0]
        )
        np.testing.assert_array_equal(
            restart.box, np.concatenate([lengths[0], angles[0]])
        )
    finally:
        restart.close()


def test_streamed_case_classifies_observed_two_leg_geometry(
    tmp_path: Path, monkeypatch: Any
) -> None:
    """Synthetic trajectory must pass both geometry legs and retain source frame."""
    from types import SimpleNamespace
    import csv

    topology = tmp_path / "topology.prmtop"
    topology.write_text("Synthetic topology ownership is explicit in test mapping.\n")
    trajectory = tmp_path / "trajectory.nc"
    xyz = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, -1, 0],
            [0, -1, 1],
            [1, -2, 0],
            [4.8, 1, 0],
            [4.8, 2, 0],
            [5.6, 0, 0],
            [6.6, 0.8, 0.4],
            [2.8, 0, 0],
            [3.8, 0, 0],
            [2.8, 1, 0],
        ],
        dtype=np.float32,
    )
    with netcdf_file(str(trajectory), "w") as data:
        for dimension, count in [("frame", 1), ("atom", len(xyz)), ("xyz", 3)]:
            data.createDimension(dimension, count)
        data.createVariable("coordinates", "f", ("frame", "atom", "xyz"))[:] = xyz[None]
        data.createVariable("cell_lengths", "d", ("frame", "xyz"))[:] = [[20, 20, 20]]
        data.createVariable("cell_angles", "d", ("frame", "xyz"))[:] = [[90, 90, 90]]
        data.createVariable("time", "d", ("frame",))[:] = [10]
    geometry = tmp_path / "geometry.csv"
    record = {
        "frame": 1,
        "time_ps": 10,
        "tier": "A",
        "projected_phi": 90,
        "d1": 3.0,
        "d2": 3.0,
        "d3": 3.0,
        "attack": 110,
    }
    with geometry.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record))
        writer.writeheader()
        writer.writerow(record)
    mapped = {
        "nalpha": 0,
        "nalpha_h": [1, 2],
        "og1": 3,
        "cb": 4,
        "carbon": 6,
        "oxygen": 7,
        "nleave": 8,
        "leaving_c1": 9,
        "waters": [(3, 10, 11, 12)],
        "ligand_bonds": [(6, 8)],
    }
    structure = SimpleNamespace(atoms=[SimpleNamespace(atomic_number=6) for _ in xyz])
    monkeypatch.setattr(SCREEN.parmed, "load_file", lambda _: structure)
    monkeypatch.setattr(SCREEN, "mapping", lambda *_: mapped)
    monkeypatch.setattr(SCREEN, "_cif_atoms", lambda *_: {})
    monkeypatch.setattr(
        SCREEN,
        "native_geometry",
        lambda *_: {
            "native_pass": True,
            "native_thr_rmsd_A": 0.5,
            "native_og1_shift_A": 0.5,
            "native_nalpha_shift_A": 0.7,
            "psi_delta_deg": 1.0,
            "chi1_delta_deg": 2.0,
        },
    )
    case = {
        "name": "fixture",
        "topology": str(topology),
        "trajectories": [str(trajectory)],
        "expected_frames": 1,
        "geometry_csv": str(geometry),
        "residue_names": {},
    }
    summary, rows = SCREEN.scan_case(case, {"reference_cif": "synthetic"}, tmp_path)
    assert summary["counts"]["joint_eligible"] == 1
    assert summary["counts"]["strict_both_leg_hbond_events"] == 1
    assert rows[0]["frame"] == 1 and rows[0]["water_resid"] == 3
    assert rows[0]["water_h_to_nleave_A"] == pytest.approx(1.8, abs=1e-6)
    assert rows[0]["relay_class"] == "strict_both_leg_hbond"
    diagnostic_path = SCRIPT.parent / "diagnose_nalpha_relay_screen.py"
    diagnostic_spec = importlib.util.spec_from_file_location(
        "relay_diagnostic", diagnostic_path
    )
    assert diagnostic_spec is not None and diagnostic_spec.loader is not None
    diagnostic = importlib.util.module_from_spec(diagnostic_spec)
    diagnostic_spec.loader.exec_module(diagnostic)
    monkeypatch.setattr(diagnostic, "mapping", lambda *_: mapped)
    diagnostic_config = tmp_path / "diagnostic_config.json"
    import json
    import sys

    diagnostic_config.write_text(json.dumps({"cases": [case]}))
    diagnostic_output = tmp_path / "diagnostic_results"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(diagnostic_path),
            "--config",
            str(diagnostic_config),
            "--primary",
            str(tmp_path),
            "--controls",
            str(tmp_path),
            "--output",
            str(diagnostic_output),
        ],
    )
    diagnostic.main()
    result = json.loads((diagnostic_output / "diagnostics.json").read_text())
    assert len(result["primary_nearpasses"]) == 1
    assert result["primary_nearpasses"][0]["max_bridge_leg_A"] == pytest.approx(
        2.8, abs=1e-6
    )
    assert result["primary"]["fixture"]["native_pass_frames"] == 1
    with pytest.raises(ValueError, match="frame count"):
        SCREEN.scan_case(
            {**case, "expected_frames": 2}, {"reference_cif": "synthetic"}, tmp_path
        )


def test_main_preserves_empty_and_separated_selection_results(
    tmp_path: Path, monkeypatch: Any
) -> None:
    """CLI must retain a negative result and select distinct frames only when present."""
    import json
    import sys

    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "source_files": [],
                "hash_inputs": [],
                "cases": [{"name": "fixture"}],
                "native_gates": {},
            }
        )
    )
    common = {
        "case": "fixture",
        "tier": "A",
        "relay_class": "leaving_leg_hbond",
        "water_h_to_nleave_A": 2.0,
        "water_h_to_nleave_angle_deg": 160.0,
        "native_thr_rmsd_A": 0.5,
        "water_resid": 3,
    }
    candidates = [{**common, "frame": frame} for frame in (1, 2, 150)]
    monkeypatch.setattr(
        SCREEN, "scan_case", lambda *_: ({"case": "fixture", "counts": {}}, candidates)
    )
    monkeypatch.setattr(SCREEN, "extract", lambda row, *_: row)
    output = tmp_path / "positive"
    monkeypatch.setattr(
        sys,
        "argv",
        [str(SCRIPT), "--config", str(config_path), "--output", str(output)],
    )
    SCREEN.main()
    assert json.loads((output / "screen.json").read_text())["selected_count"] == 2
    monkeypatch.setattr(
        SCREEN, "scan_case", lambda *_: ({"case": "fixture", "counts": {}}, [])
    )
    output = tmp_path / "negative"
    monkeypatch.setattr(
        sys,
        "argv",
        [str(SCRIPT), "--config", str(config_path), "--output", str(output)],
    )
    SCREEN.main()
    assert json.loads((output / "screen.json").read_text())["selected_count"] == 0
