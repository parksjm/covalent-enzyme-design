"""Focused regressions for the approved63atom donor cluster construction."""

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from rdkit import Chem

import prepare_oxyanion_cluster as builder


def test_inventory_roundtrip_and_provenance(tmp_path: Path, monkeypatch: Any) -> None:
    """Preserve all42 input coordinates in both outputs and retain hashed sources."""
    shutil.copy2(builder.BASE / "donor_fragments.py", tmp_path / "donor_fragments.py")
    monkeypatch.setattr(builder, "BASE", tmp_path)
    builder.main()
    out = tmp_path / "seed_01"
    molecule = Chem.SDMolSupplier(str(out / "native_seed.sdf"), removeHs=False)[0]
    assert molecule is not None
    xyz = np.asarray(
        [
            [float(v) for v in s.split()[1:]]
            for s in (out / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    np.testing.assert_array_equal(molecule.GetConformer().GetPositions(), xyz)
    np.testing.assert_array_equal(xyz[:42], builder._read_seed()[1])
    records = json.loads((out / "coordinate_provenance.json").read_text())
    assert len(records) == 63
    assert sum(row["constructed_cap_H"] for row in records[42:]) == 5
    report = json.loads((out / "preflight.json").read_text())
    assert report["formula"] == "C15H37N5O6"
    assert report["component_atom_counts"] == [3, 9, 12, 39]
    assert report["geometry_screen_passed"]
    assert report["coordinate_stereocenters"] == {
        "THR:CA": "S",
        "THR:CB": "R",
        "PA6:C54": "R",
    }
    for row in json.loads((out / "source_manifest.json").read_text()):
        assert (
            hashlib.sha256((out / row["snapshot"]).read_bytes()).hexdigest()
            == row["sha256"]
        )
    with pytest.raises(FileExistsError):
        builder.main()


def test_contact_and_stereo_corruption_are_rejected() -> None:
    """Reject a donor proton collision and inversion of the native stereochemistry."""
    molecule, xyz, _, _ = builder.build_cluster()
    indices = {a.GetProp("origin_id"): a.GetIdx() for a in molecule.GetAtoms()}
    collision = xyz.copy()
    collision[indices["NMA:113:H"]] = xyz[indices["PA6:O49"]] + [0.1, 0.0, 0.0]
    report = builder.assess_cluster(molecule, collision)
    assert not report["geometry_screen_passed"]
    assert report["geometry"]["geometry_failures"]
    mirrored = xyz.copy()
    mirrored[:, 0] *= -1
    conformer = molecule.GetConformer()
    for index, point in enumerate(mirrored):
        conformer.SetAtomPosition(index, tuple(float(v) for v in point))
    with pytest.raises(ValueError, match="Coordinate-derived stereochemistry mismatch"):
        builder.assess_cluster(molecule, mirrored)


def test_misaligned_source_frame_fails(monkeypatch: Any) -> None:
    """Fail explicitly if the seed and donor source cease sharing a frame."""
    original = builder._read_seed()
    monkeypatch.setattr(
        builder, "_read_seed", lambda: (original[0], original[1] + 0.1, original[2])
    )
    with pytest.raises(ValueError, match="Source coordinate-frame mismatch"):
        builder.build_cluster()
