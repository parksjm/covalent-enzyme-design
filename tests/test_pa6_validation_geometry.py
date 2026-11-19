"""Test strict coordinate-only transfer into typed MOL2 files."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.pa6_validation_geometry import (
    prepare_pa6_validation_geometry,
    replace_mol2_coordinates,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "provenance/xa_ma_capped_reference_20261022/pa6_l4_source.mol2"


MOL2 = """@<TRIPOS>MOLECULE
test
 2 1 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 A1 0.0 0.0 0.0 c3 1 LIG -0.1
2 A2 1.0 0.0 0.0 hc 1 LIG 0.1
@<TRIPOS>BOND
1 1 2 1
"""


def test_replace_mol2_coordinates_preserves_types_and_charges() -> None:
    rendered = replace_mol2_coordinates(
        MOL2, ["A1", "A2"], np.asarray([[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]])
    )

    assert "1 A1 2.0000000000 3.0000000000 4.0000000000 c3 1 LIG -0.1" in rendered
    assert "2 A2 5.0000000000 6.0000000000 7.0000000000 hc 1 LIG 0.1" in rendered


@pytest.mark.parametrize("names", [["A2", "A1"], ["A1"]])
def test_replace_mol2_coordinates_rejects_identity_mismatch(names: list[str]) -> None:
    with pytest.raises(ValueError):
        replace_mol2_coordinates(MOL2, names, np.zeros((len(names), 3)))


def test_prepare_validation_geometry_records_accepted_source(tmp_path: Path) -> None:
    molecule = load_pa6(SOURCE).mol
    xyz = np.asarray(molecule.GetConformer().GetPositions())
    rows = [
        f"{index:6d}{atom.GetAtomicNum():11d}{0:12d}"
        f"{point[0]:16.6f}{point[1]:12.6f}{point[2]:12.6f}"
        for index, (atom, point) in enumerate(
            zip(molecule.GetAtoms(), xyz, strict=True), start=1
        )
    ]
    gaussian = tmp_path / "optimization.log"
    gaussian.write_text(
        "\n".join(
            [
                "Standard orientation:",
                "---------------------------------------------------------------------",
                "Center Atomic Atomic Coordinates (Angstroms)",
                "Number Number Type X Y Z",
                "---------------------------------------------------------------------",
                *rows,
                "---------------------------------------------------------------------",
                "Normal termination of Gaussian 16",
                "",
            ]
        )
    )
    acceptance = tmp_path / "acceptance.json"
    acceptance.write_text(json.dumps({"accepted": True, "conformer": "test"}))
    output = tmp_path / "typed_accepted.mol2"
    manifest = tmp_path / "manifest.json"

    report = prepare_pa6_validation_geometry(
        SOURCE, SOURCE, gaussian, acceptance, output, manifest
    )

    assert report["atom_count"] == 79
    assert report["conformer"] == "test"
    assert output.is_file()
    assert json.loads(manifest.read_text()) == report


def test_prepare_validation_geometry_requires_acceptance(tmp_path: Path) -> None:
    acceptance = tmp_path / "acceptance.json"
    acceptance.write_text(json.dumps({"accepted": False}))
    with pytest.raises(RuntimeError, match="positive scientific acceptance"):
        prepare_pa6_validation_geometry(
            SOURCE,
            SOURCE,
            tmp_path / "missing.log",
            acceptance,
            tmp_path / "output.mol2",
            tmp_path / "manifest.json",
        )
