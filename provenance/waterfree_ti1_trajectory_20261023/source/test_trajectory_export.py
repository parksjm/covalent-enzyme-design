"""Test complete and partial live-log trajectory exports without fabricated frames."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from rdkit import Chem

BASE = Path(__file__).resolve().parent
PILOT = BASE / "waterfree_pilot_01"


def module() -> ModuleType:
    """Load the retained visualization exporter."""
    spec = importlib.util.spec_from_file_location(
        "trajectory_export", BASE / "export_optimization_trajectory.py"
    )
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def fixture_log(
    path: Path,
    *,
    broken: bool = False,
    partial: bool = False,
    wrong_element: bool = False,
) -> bytes:
    """Write explicitly synthetic Gaussian tables for physical corruption tests."""
    mol = Chem.SDMolSupplier(str(PILOT / "native_seed.sdf"), removeHs=False)[0]
    xyz = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (PILOT / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    mapping = json.loads((PILOT / "atom_map.json").read_text())["atom_index_by_origin"]
    lines = ["SYNTHETIC TEST FIXTURE; NOT AN ENGINE RESULT"]
    for frame in range(2):
        points = xyz.copy()
        if frame:
            points[mapping["THR:H1"], 0] += 0.01
            if broken:
                points[mapping["PA6:C54"], 0] += 4
        lines += [
            " Input orientation:",
            " -----",
            " Center Atomic Atomic Coordinates (Angstroms)",
            " Number Number Type X Y Z",
            " -----",
        ]
        for i, (atom, point) in enumerate(zip(mol.GetAtoms(), points, strict=True), 1):
            number = 6 if wrong_element and i == 1 else atom.GetAtomicNum()
            lines += [f"{i} {number} 0 " + " ".join(f"{v:.6f}" for v in point)]
        lines += [" -----"]
    if partial:
        lines += [
            " Input orientation:",
            " -----",
            " Center Atomic Atomic Coordinates (Angstroms)",
        ]
    data = ("\n".join(lines) + "\n").encode()
    path.write_bytes(data)
    return data


def invoke(helper: ModuleType, log: Path, output: Path) -> dict:
    """Use the public exporter with the exact matching chemical graph."""
    return helper.export(
        log,
        PILOT / "native_seed.sdf",
        PILOT / "atom_map.json",
        PILOT / "pilot_policy.json",
        output,
    )


def test_multi_frame_export_preserves_hydrogens(tmp_path: Path) -> None:
    """Require real table count, exact snapshot bytes and portable movie files."""
    log = tmp_path / "synthetic.log"
    data = fixture_log(log)
    output = tmp_path / "export"
    report = invoke(module(), log, output)
    assert report["frames"] == 2
    assert report["atoms_per_frame"] == 39
    assert report["explicit_H_per_frame"] == 23
    assert report["source_sha256"] == hashlib.sha256(data).hexdigest()
    assert report["geometry_acceptance_claimed"] is False
    assert not report["optimization_completed"]
    assert not report["normal_termination"]
    assert (output / "optimization.log.txt").read_bytes() == data
    assert log.read_bytes() == data
    pdb = (output / "optimization_trajectory.pdb").read_text()
    assert sum(line.startswith("MODEL ") for line in pdb.splitlines()) == 2
    assert sum(line.startswith(("ATOM  ", "HETATM")) for line in pdb.splitlines()) == 78
    assert "CONECT" in pdb
    mol = Chem.MolFromPDBFile(
        str(output / "optimization_trajectory.pdb"),
        removeHs=False,
        sanitize=False,
        proximityBonding=False,
    )
    assert mol is not None and mol.GetNumAtoms() == 39
    assert sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) == 23
    source_mol = Chem.SDMolSupplier(str(PILOT / "native_seed.sdf"), removeHs=False)[0]
    assert mol.GetNumBonds() == source_mol.GetNumBonds()
    assert [a.GetFormalCharge() for a in mol.GetAtoms()] == [
        a.GetFormalCharge() for a in source_mol.GetAtoms()
    ]
    assert sum(b.GetBondTypeAsDouble() for b in mol.GetBonds()) == sum(
        b.GetBondTypeAsDouble() for b in source_mol.GetBonds()
    )
    viewer = (output / "view_optimization.pml").read_text()
    assert "set valence, 0" in viewer
    assert "set movie_fps, 3" in viewer
    assert "mset 1 -2" in viewer
    assert not any(line.strip() == "python" for line in viewer.splitlines())
    assert (output / "trajectory_metrics.csv").is_file()
    assert (output / "trajectory_report.json").is_file()


def test_partial_table_not_invented(tmp_path: Path) -> None:
    """Export complete frames while explicitly reporting an incomplete last table."""
    log = tmp_path / "partial.log"
    fixture_log(log, partial=True)
    report = invoke(module(), log, tmp_path / "export")
    assert report["frames"] == 2
    assert report["ignored_trailing_incomplete_counts"]["Input"] == 1


def test_bad_geometry_still_available_for_inspection(tmp_path: Path) -> None:
    """Do not hide a broken TI1 structure by dropping its trajectory frame."""
    log = tmp_path / "broken.log"
    fixture_log(log, broken=True)
    output = tmp_path / "export"
    report = invoke(module(), log, output)
    assert report["frames"] == 2
    assert (output / "optimization_trajectory.xyz").read_text().count("39\n") == 2
    assert report["geometry_acceptance_claimed"] is False


def test_wrong_elements_fail(tmp_path: Path) -> None:
    """Fail instead of mapping a log to a mismatched chemical topology."""
    log = tmp_path / "wrong.log"
    fixture_log(log, wrong_element=True)
    with pytest.raises(ValueError):
        invoke(module(), log, tmp_path / "export")


def test_refuses_overwrite(tmp_path: Path) -> None:
    """An export refresh cannot overwrite its source or an earlier snapshot."""
    log = tmp_path / "source.log"
    fixture_log(log)
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("keep")
    with pytest.raises(FileExistsError):
        invoke(module(), log, output)
    assert marker.read_text() == "keep"
