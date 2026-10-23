"""Verify truncation preserves native chemistry and exact retained coordinates."""

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
    """Load the retained calculation preparer without running its entry point."""
    spec = importlib.util.spec_from_file_location(
        "truncated_preparer", BASE / "prepare_truncated_pilot.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_real_source_preserves_local_chemistry(tmp_path: Path) -> None:
    """Check cuts, proton inventory, state, mapping and Gaussian resource inputs."""
    module = load_preparer()
    module.OUT = tmp_path / "pilot"
    module.main()
    molecule = Chem.SDMolSupplier(str(module.OUT / "native_seed.sdf"), removeHs=False)[
        0
    ]
    metadata = json.loads((module.OUT / "atom_map.json").read_text())
    mapping = metadata["atom_index_by_origin"]
    summary = json.loads((module.OUT / "preflight.json").read_text())
    assert summary["accepted"]
    assert summary["atom_count"] == 42
    assert summary["explicit_H"] == 25
    assert summary["retained_coordinate_max_displacement_A"] == 0
    assert Chem.GetFormalCharge(molecule) == 0
    assert all(a.GetNumImplicitHs() == 0 for a in molecule.GetAtoms())
    for center in ("PA6:C70", "PA6:C75"):
        atom = molecule.GetAtomWithIdx(mapping[center])
        assert sum(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()) == 3
    for center in ("PA6:C62", "PA6:C67"):
        atom = molecule.GetAtomWithIdx(mapping[center])
        assert sum(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()) == 2
    thr_n = molecule.GetAtomWithIdx(mapping["THR:N"])
    assert thr_n.GetFormalCharge() == 1
    assert sum(n.GetAtomicNum() == 1 for n in thr_n.GetNeighbors()) == 3
    assert molecule.GetAtomWithIdx(mapping["PA6:O49"]).GetFormalCharge() == -1
    assert molecule.GetAtomWithIdx(mapping["PA6:N59"]).GetFormalCharge() == 0
    assert molecule.GetAtomWithIdx(mapping["PA6:C54"]).GetProp("_CIPCode") == "R"
    assert molecule.GetAtomWithIdx(mapping["WAT:O"]).GetDegree() == 2
    assert "PA6:N56" not in mapping
    assert "PA6:C53" not in mapping
    xyz = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (module.OUT / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    for center, cap in [("PA6:C70", "PA6:HC70"), ("PA6:C75", "PA6:HC75")]:
        assert np.linalg.norm(
            xyz[mapping[center]] - xyz[mapping[cap]]
        ) == pytest.approx(1.09, abs=1e-9)
    policy = json.loads((module.OUT / "pilot_policy.json").read_text())
    assert len(policy["frozen_atom_indices_0based"]) == 12
    assert policy["nprocshared"] == 16
    text = (module.OUT / "input.gjf").read_text()
    assert "%nprocshared=16" in text
    assert "0 1" in text
    assert "Freq" not in text
    assert "CalcFC" not in text


def test_rejects_wrong_source_inventory(tmp_path: Path) -> None:
    """Prevent accidental recursive truncation of an already truncated model."""
    module = load_preparer()
    shutil.copy2(
        BASE / "truncated_pilot_01/native_seed.sdf", tmp_path / "native_seed.sdf"
    )
    module.SOURCE = tmp_path
    module.OUT = tmp_path / "pilot"
    with pytest.raises(ValueError, match="103-atom"):
        module.main()
    assert not module.OUT.exists()
