"""Deliberately corrupt native coordinates/state to test the pilot gate."""

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
SOURCE = Path(
    "/data/jpc/covalent-enzyme-design/output_nyl12_smd_water_bridge_20261023/svp_pilot_01"
)


def load_checker() -> ModuleType:
    """Load the retained assessment script for isolated fixture testing."""
    spec = importlib.util.spec_from_file_location(
        "native_pilot_assessor", SOURCE / "assess_pilot.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    "fault",
    [
        "none",
        "broken_reactive_bond",
        "frozen_atom_moved",
        "wrong_charge",
        "not_converged",
        "proton_transfer",
        "proton_return_OG1",
        "water_bridge_lost",
    ],
)
def test_rejects_invalid_native_states(tmp_path: Path, fault: str) -> None:
    """Reject physical/state corruptions rather than trusting normal exit text."""
    source = SOURCE
    for name in (
        "native_seed.sdf",
        "native_seed.xyz",
        "atom_map.json",
        "pilot_policy.json",
    ):
        shutil.copy2(source / name, tmp_path / name)
    xyz = np.asarray(
        [
            [float(v) for v in line.split()[1:]]
            for line in (source / "native_seed.xyz").read_text().splitlines()[2:]
        ]
    )
    origins = json.loads((source / "atom_map.json").read_text())["atom_index_by_origin"]
    if fault == "broken_reactive_bond":
        xyz[origins["PA6:C54"]] += [4, 0, 0]
    elif fault == "frozen_atom_moved":
        xyz[origins["THR:CA"]] += [0.2, 0, 0]
    elif fault == "proton_transfer":
        xyz[origins["THR:H1"]] = xyz[origins["PA6:O49"]] + [0.95, 0, 0]
    elif fault == "proton_return_OG1":
        xyz[origins["THR:HG1"]] = xyz[origins["THR:OG1"]] + [0.95, 0, 0]
    elif fault == "water_bridge_lost":
        for name in ("WAT:O", "WAT:H1", "WAT:H2"):
            xyz[origins[name]] += [0, 0, 4]
    molecule = Chem.SDMolSupplier(str(source / "native_seed.sdf"), removeHs=False)[0]
    lines = [
        "SYNTHETIC UNIT TEST; NOT A GAUSSIAN RESULT",
        "Charge = 1 Multiplicity = 1"
        if fault == "wrong_charge"
        else "Charge = 0 Multiplicity = 1",
        "Input orientation:",
        "-----",
        "headers",
        "-----",
    ]
    for idx, (atom, point) in enumerate(zip(molecule.GetAtoms(), xyz, strict=True), 1):
        lines.append(
            f"{idx} {atom.GetAtomicNum()} 0 " + " ".join(f"{v:.10f}" for v in point)
        )
    lines += ["-----"]
    if fault != "not_converged":
        lines += ["Optimization completed."]
    lines += ["Normal termination of Gaussian 16 synthetic test fixture"]
    (tmp_path / "optimization.log").write_text("\n".join(lines) + "\n")
    checker = load_checker()
    checker.__file__ = str(tmp_path / "assess_pilot.py")
    try:
        checker.main()
    except SystemExit as exc:
        assert exc.code == 3
    report = json.loads((tmp_path / "scientific_acceptance.json").read_text())
    assert report["accepted"] == all(report["checks"].values())
    assert report["forcefield_accepted"] is False
    assert report["charge_fit_reference_accepted"] is False
    assert (tmp_path / "final_structure.xyz").is_file()
    if fault == "broken_reactive_bond":
        assert not report["checks"]["reaction_center"]
        assert not report["checks"]["graph_and_stereo"]
    elif fault == "frozen_atom_moved":
        assert not report["checks"]["fixed_coordinates"]
    elif fault == "wrong_charge":
        assert not report["checks"]["electronic_state"]
    elif fault == "not_converged":
        assert not report["checks"]["optimized"]
    elif fault == "proton_transfer":
        assert not report["checks"]["oxyanion_not_protonated_screen"]
    elif fault == "proton_return_OG1":
        assert not report["checks"]["ammonium_retained"]
    elif fault == "water_bridge_lost":
        assert not report["checks"]["water_bridge"]
    else:
        assert report["checks"]["water_bridge"]
        assert report["checks"]["ammonium_retained"]
        assert report["accepted"]
        assert report["checks"]["fixed_coordinates"]
        assert report["checks"]["chi1"]
        assert report["checks"]["electronic_state"]
        assert report["checks"]["optimized"]
    if fault != "none":
        assert not report["accepted"]
