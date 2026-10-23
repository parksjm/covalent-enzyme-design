"""Test water-only removal and unchanged physical constraints in GIC input."""

from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from rdkit import Chem

BASE = Path(__file__).resolve().parent


def load_preparer() -> ModuleType:
    """Load the retained water-free preparation module without running main."""
    spec = importlib.util.spec_from_file_location(
        "waterfree_preparer", BASE / "prepare_waterfree39.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_water_only_removal_preserves_chemistry(tmp_path: Path) -> None:
    """Preserve all39 atoms, charges, caps, stereo and physical freeze choices."""
    module = load_preparer()
    module.OUT = tmp_path / "candidate"
    module.main()
    mol = Chem.SDMolSupplier(str(module.OUT / "native_seed.sdf"), removeHs=False)[0]
    metadata = json.loads((module.OUT / "atom_map.json").read_text())
    mapping = metadata["atom_index_by_origin"]
    assert mol.GetNumAtoms() == 39
    assert Chem.GetFormalCharge(mol) == 0
    assert sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) == 23
    assert all(a.GetNumImplicitHs() == 0 for a in mol.GetAtoms())
    assert not any(origin.startswith("WAT:") for origin in mapping)
    assert mol.GetAtomWithIdx(mapping["THR:N"]).GetFormalCharge() == 1
    assert (
        sum(
            n.GetAtomicNum() == 1
            for n in mol.GetAtomWithIdx(mapping["THR:N"]).GetNeighbors()
        )
        == 3
    )
    assert mol.GetAtomWithIdx(mapping["PA6:O49"]).GetFormalCharge() == -1
    assert mol.GetAtomWithIdx(mapping["PA6:C54"]).GetProp("_CIPCode") == "R"
    source = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (module.SOURCE / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    target = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (module.OUT / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    provenance = json.loads((module.OUT / "coordinate_provenance.json").read_text())
    assert np.array_equal(
        target, source[[row["source_42atom_index_0based"] for row in provenance]]
    )
    oldmap = json.loads((module.SOURCE / "atom_map.json").read_text())[
        "atom_index_by_origin"
    ]
    old_inverse = {i: origin for origin, i in oldmap.items()}
    oldpolicy = json.loads((module.SOURCE / "pilot_policy.json").read_text())
    policy = json.loads((module.OUT / "pilot_policy.json").read_text())
    assert {old_inverse[i] for i in oldpolicy["frozen_atom_indices_0based"]} == {
        origin
        for origin, i in mapping.items()
        if i in policy["frozen_atom_indices_0based"]
    }
    text = (module.OUT / "input.gjf").read_text()
    assert text.count("(Frozen)=") == 37
    assert text.count("=X(") == text.count("=Y(") == text.count("=Z(") == 39
    assert "CHI(Frozen)=D(1,2,3,4)" in text
    assert "ReadAllGIC" in text and "MaxStep=2" in text
    assert "ModRedundant" not in text
    assert "%nprocshared=16" in text
    assert "Freq" not in text and "CalcFC" not in text
    assert "WAT:O" not in policy["mobile_heavy_origins"]


def test_rejects_wrong_source_inventory(tmp_path: Path) -> None:
    """Reject recursive removal from an already water-free39 atom model."""
    module = load_preparer()
    shutil.copy2(
        BASE / "waterfree_pilot_01/native_seed.sdf", tmp_path / "native_seed.sdf"
    )
    module.SOURCE = tmp_path
    module.OUT = tmp_path / "candidate"
    with pytest.raises(ValueError, match="42-atom"):
        module.main()
    assert not module.OUT.exists()
