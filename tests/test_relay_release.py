"""Scientific input and geometry contracts for native relay release jobs."""

from pathlib import Path

import numpy as np
import pytest
from nyl12_covalent.relay_release import (
    ROLES,
    aligned_rmsd,
    geometry,
    minimization_input,
)


def test_rmsd_removes_rigid_motion_preserves_distortion() -> None:
    xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    rotation = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert aligned_rmsd(xyz, xyz @ rotation + 14) < 1e-14
    changed = xyz.copy()
    changed[0, 0] += 0.3
    assert aligned_rmsd(xyz, changed) > 0.05


@pytest.mark.parametrize(
    "other", [np.zeros((3, 2)), np.empty((0, 3)), np.full((3, 3), np.nan)]
)
def test_rmsd_rejects_invalid(other: np.ndarray) -> None:
    with pytest.raises(ValueError):
        aligned_rmsd(other, other)


def test_geometry_records_all_proton_alternatives() -> None:
    xyz = np.ones((max(ROLES.values()), 3)) * 20
    positions = {
        "C": [0, 0, 0],
        "O": [0, 0, 1.3],
        "OG": [1.5, 0, 0],
        "NL": [0, 1.5, 0],
        "CA": [-1.5, 0, 0],
        "N": [5, 5, 5],
        "HG": [5, 5, 6],
        "H1": [5, 6, 5],
        "H2": [6, 5, 5],
        "HN": [0, 2.5, 0],
    }
    for role, value in positions.items():
        xyz[ROLES[role] - 1] = value
    result = geometry(xyz)
    assert result["geometry_label"] == "TI1_like_geometric_candidate"
    assert result["proton_bonds"] == {
        "HG": ["N"],
        "H1": ["N"],
        "H2": ["N"],
        "HN": ["NL"],
    }
    assert len(result["carbon_angles_degrees"]) == 6
    xyz[ROLES["OG"] - 1] = [3.0, 0, 0]
    assert geometry(xyz)["geometry_label"] == "reactant_like"
    xyz[ROLES["OG"] - 1] = [1.5, 0, 0]
    xyz[ROLES["NL"] - 1] = [0, 3, 0]
    assert geometry(xyz)["geometry_label"] == "acyl_like"
    xyz[ROLES["NL"] - 1] = [0, 1.8, 0]
    assert geometry(xyz)["geometry_label"] == "unassigned"
    xyz[ROLES["OG"] - 1] = 0
    with pytest.raises(ValueError, match="Coincident"):
        geometry(xyz)


@pytest.mark.parametrize("xyz", [np.zeros((2, 3)), np.full((19937, 3), np.nan)])
def test_geometry_rejects_invalid(xyz: np.ndarray) -> None:
    with pytest.raises(ValueError):
        geometry(xyz)


def test_native_input_preserves_boundary_and_hamiltonian() -> None:
    source = """&cntrl\n ntr=0, nmropt=0, maxcyc=10000,\n/\n&qmmm\n dftb_telec=0.0, qmmask='@1-4', qmcharge=-2, qm_theory='DFTB3',\n qm_ewald=1, qm_pme=1, qmcut=10.0,\n/\n&dlfind\n active='@1-7', dimer=.true.,\n/\n"""
    result = minimization_input(source)
    for text in (
        "active='@1-7'",
        "qmmask='@1-4'",
        "qmcharge=-2",
        "dimer=.false.",
        "maxcycle=10000",
        "tol=1.d-5",
        "ntr=0",
    ):
        assert text in result
    assert "restraint" not in result
    with pytest.raises(ValueError):
        minimization_input(source, 0)
    with pytest.raises(ValueError):
        minimization_input(source.replace("ntr=0", "ntr=1"))
    with pytest.raises(ValueError):
        minimization_input(source.replace("nmropt=0", "nmropt=1"))


@pytest.mark.parametrize(
    "converged,exit_code,force",
    [(True, 0, 0.0), (False, 0, 0.0), (True, 1, 0.0), (True, 0, 1.0)],
)
def test_native_execution_acceptance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    converged: bool,
    exit_code: int,
    force: float,
) -> None:
    import json
    from types import SimpleNamespace
    from nyl12_covalent import relay_release_run as module

    seed = tmp_path / "seed"
    seed.mkdir()
    (tmp_path / "prepared.json").write_text(json.dumps({"parameter_sha256": {}}))
    xyz = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    class State:
        def __init__(self, root: Path) -> None:
            self.xyz = xyz
            self.sander = Path("/fake/sander")
            self.active = np.array([True, False])
            self.count = 0

        def gradient(
            self, coordinates: np.ndarray, evidence: Path, *, precision: bool
        ) -> tuple[float, np.ndarray]:
            self.count += 1
            return -float(self.count), np.zeros((1, 3)) + force

    def process(command: list[str], **kwargs: object) -> SimpleNamespace:
        (seed / "minimize.mdout").write_text(
            "Optimisation converged" if converged else "limit"
        )
        return SimpleNamespace(returncode=exit_code)

    monkeypatch.setattr(module, "NativeState", State)
    monkeypatch.setattr(module.subprocess, "run", process)
    monkeypatch.setattr(module, "audit_dftb_scc", lambda text: {"passed": True})
    monkeypatch.setattr(module, "geometry", lambda x: {"test_geometry": True})
    monkeypatch.setattr(
        module.Rst7, "open", lambda path: SimpleNamespace(coordinates=xyz)
    )
    if converged and exit_code == 0 and force == 0:
        module.run(seed)
        result = json.loads((seed / "scientific_acceptance.json").read_text())
        assert result["accepted"] and result["frozen_max_displacement_A"] == 0
        assert result["energy_decreased"]
    else:
        with pytest.raises(RuntimeError):
            module.run(seed)
        assert not json.loads((seed / "scientific_acceptance.json").read_text())[
            "accepted"
        ]
    assert (
        json.loads((seed / "process_acceptance.json").read_text())["exit_status"]
        == exit_code
    )
