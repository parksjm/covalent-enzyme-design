"""Test strict generation of bidirectional PA6 relaxed scans."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from rdkit import Chem

from nyl12_covalent.pa6_torsion_fragments import prepare_torsion_fragments
from nyl12_covalent.pa6_torsion_scans import SCANS, prepare_relaxed_scans


def _write_completed_logs(geometry_dir: Path) -> None:
    manifest = json.loads((geometry_dir / "manifest.json").read_text())
    for record in manifest["fragments"]:
        fragment = geometry_dir / record["label"]
        molecule = next(
            item
            for item in Chem.SDMolSupplier(
                str(fragment / "seed_mmff94s.sdf"), removeHs=False
            )
            if item is not None
        )
        coordinates = molecule.GetConformer().GetPositions()
        rows = [
            f"{index:6d}{atom.GetAtomicNum():11d}{0:12d}"
            f"{point[0]:16.6f}{point[1]:12.6f}{point[2]:12.6f}"
            for index, (atom, point) in enumerate(
                zip(molecule.GetAtoms(), coordinates, strict=True), start=1
            )
        ]
        (fragment / "optimization.log").write_text(
            "\n".join(
                [
                    "Standard orientation:",
                    "---------------------------------------------------------------------",
                    "Center Atomic Atomic Coordinates (Angstroms)",
                    "Number Number Type X Y Z",
                    "---------------------------------------------------------------------",
                    *rows,
                    "---------------------------------------------------------------------",
                    "Optimization completed.",
                    "Normal termination of Gaussian 16",
                    "",
                ]
            )
        )


def test_prepare_bidirectional_scans(tmp_path: Path) -> None:
    geometries = tmp_path / "geometries"
    prepare_torsion_fragments(geometries)
    _write_completed_logs(geometries)
    output = tmp_path / "scans"

    manifest = prepare_relaxed_scans(geometries, output)

    assert manifest["scan_count"] == len(SCANS) == 8
    assert manifest["gaussian_job_count"] == 16
    for scan in manifest["scans"]:
        for direction, increment in (("forward", 1.0), ("reverse", -1.0)):
            path = Path(scan["directions"][direction]["input"])
            text = path.read_text()
            assert "Opt=(ModRedundant,Tight,MaxCycles=200)" in text
            assert "CalcFC" not in text
            assert (
                f" {increment * abs(scan['directions'][direction]['increment_degrees']):.6f}"
                in text
            )


def test_scan_preparation_requires_completed_geometry(tmp_path: Path) -> None:
    geometries = tmp_path / "geometries"
    prepare_torsion_fragments(geometries)
    with pytest.raises(RuntimeError, match="optimization incomplete"):
        prepare_relaxed_scans(geometries, tmp_path / "scans")
