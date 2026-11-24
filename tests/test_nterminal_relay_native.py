"""Mocked native execution and scientific acceptance contracts."""

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from nyl12_covalent import nterminal_relay_native as native


@pytest.mark.parametrize(
    ("text", "biased", "physical"),
    [
        ("Etot = -20.0 RESTRAINT = 3.0", True, -23.0),
        ("Etot = 4D+00", False, 4.0),
        ("Etot = 1 Etot = 2", False, 2.0),
    ],
)
def test_energies(text: str, biased: bool, physical: float) -> None:
    assert native.native_energies(text, biased)["physical_kcal_mol"] == physical


@pytest.mark.parametrize(
    "text,biased",
    [
        ("none", False),
        ("Etot = 1", True),
        ("Etot = 1 RESTRAINT = -2", True),
        ("Etot = 1e999", False),
    ],
)
def test_invalid_energies(text: str, biased: bool) -> None:
    with pytest.raises(ValueError):
        native.native_energies(text, biased)


def test_exterior() -> None:
    xyz = np.zeros((2, 3))
    active = np.array([True, False])
    final = xyz.copy()
    final[0, 0] = 2
    assert native.unchanged_exterior(xyz, final, active) == 0
    final[1, 0] = 1
    with pytest.raises(ValueError, match="exterior"):
        native.unchanged_exterior(xyz, final, active)
    with pytest.raises(ValueError):
        native.unchanged_exterior(xyz, xyz[:1], active)
    final[1, 0] = np.nan
    with pytest.raises(ValueError):
        native.unchanged_exterior(xyz, final, active)


def test_hydrogen_audit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    model = tmp_path / "model"
    (model / "inputs").mkdir(parents=True)
    (model / "inputs/dimer.mdin").write_text("qmmask='@1'\n")
    parent = SimpleNamespace(idx=0, atomic_number=6, bond_partners=[])
    hydrogen = SimpleNamespace(idx=1, atomic_number=1, bond_partners=[parent])
    excluded = SimpleNamespace(idx=13905, atomic_number=1, bond_partners=[])
    monkeypatch.setattr(
        native,
        "AmberParm",
        lambda path: SimpleNamespace(atoms=[parent, hydrogen, excluded]),
    )
    monkeypatch.setattr(
        native,
        "AmberMask",
        lambda parm, mask: SimpleNamespace(Selection=lambda: [1, 1, 1]),
    )
    xyz = np.zeros((13906, 3))
    xyz[1, 0] = 1.0
    assert native.hydrogen_audit(model, xyz)["passed"]
    xyz[1, 0] = 2.0
    assert not native.hydrogen_audit(model, xyz)["passed"]
    hydrogen.bond_partners = []
    with pytest.raises(ValueError, match="parent"):
        native.hydrogen_audit(model, xyz)
    (model / "inputs/dimer.mdin").write_text("no mask")
    with pytest.raises(ValueError, match="QM mask"):
        native.hydrogen_audit(model, xyz)


MDIN = """&cntrl ntr=0, nmropt=0, /
&qmmm qmmask='@1', qm_theory='DFTB3', qmcharge=-2, scfconv=1.0d-12,
 dftb_telec=0.0, qmcut=10.0, qm_ewald=1, qm_pme=1, /
&dlfind active='@1', dimer=.false., /
"""


def setup_native(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, Mock, Mock, dict]:
    """Create a fake native environment with configurable engine behavior."""
    campaign = tmp_path / "campaign"
    model = campaign / "model"
    (model / "inputs").mkdir(parents=True)
    (model / "library_build").mkdir()
    (model / "library_build/libdlfind.so").write_text("lib")
    topology = model / "inputs/system.prmtop"
    topology.write_text("top")
    start = model / "start.ncrst"
    start.write_text("start")
    engine = tmp_path / "sander"
    engine.write_text("engine")
    (campaign / "config.json").write_text(
        json.dumps(
            {"model": "model", "source_topology_sha256": native.digest(topology)}
        )
    )
    state = SimpleNamespace(
        topology=topology,
        mdin_source=MDIN,
        active=np.array([True, False]),
        xyz=np.zeros((2, 3)),
        sander=engine,
        gradient=Mock(return_value=(-1.0, np.array([[1e-6, 0.0, 0.0]]))),
    )
    monkeypatch.setattr(native, "NativeState", lambda model: state)
    monkeypatch.setattr(
        native.Rst7, "open", lambda path: SimpleNamespace(coordinates=np.zeros((2, 3)))
    )
    monkeypatch.setattr(native, "action", Mock())
    monkeypatch.setattr(native, "geometry", lambda xyz: {"valid": True})
    monkeypatch.setattr(native, "hydrogen_audit", lambda model, xyz: {"passed": True})
    scc = Mock(return_value={"passed": True})
    monkeypatch.setattr(native, "audit_dftb_scc", scc)
    settings = {
        "exit": 0,
        "text": "Optimisation converged\nEtot = -10\nRESTRAINT = 2\n",
    }

    def engine_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        Path(command[command.index("-o") + 1]).write_text(settings["text"])
        Path(command[command.index("-r") + 1]).write_text("restart")
        return SimpleNamespace(returncode=settings["exit"])

    engine_mock = Mock(side_effect=engine_run)
    monkeypatch.setattr(native.subprocess, "run", engine_mock)
    return campaign, start, engine_mock, state.gradient, settings


def test_native_acceptance_and_independent_force(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign, start, engine, gradient, settings = setup_native(tmp_path, monkeypatch)
    result = native.run_native(campaign, campaign / "rel", start, MDIN)
    assert (
        result["accepted"]
        and result["physical_kcal_mol"] == -10
        and gradient.call_count == 1
    )
    gradient.return_value = (-1.0, np.array([[1e-3, 0.0, 0.0]]))
    assert not native.run_native(campaign, campaign / "high_force", start, MDIN)[
        "accepted"
    ]
    ts_input = MDIN.replace("dimer=.false.", "dimer=.true.")
    native.run_native(campaign, campaign / "ts", start, ts_input)
    assert gradient.call_count == 2
    biased = MDIN.replace("nmropt=0", "nmropt=1")
    result = native.run_native(
        campaign, campaign / "drive", start, biased, restraints="rst"
    )
    assert result["physical_kcal_mol"] == -12 and gradient.call_count == 2
    assert engine.call_count == 4
    settings["text"] = "Etot = 2"
    assert not native.run_native(campaign, campaign / "nonconverged", start, MDIN)[
        "accepted"
    ]


@pytest.mark.parametrize(
    "old,new",
    [
        ("qmmask='@1'", "qmmask='@2'"),
        ("active='@1'", "active='@2'"),
        ("qmcharge=-2", "qmcharge=-1"),
        ("scfconv=1.0d-12", "scfconv=1.0d-8"),
        ("qm_pme=1", "qm_pme=0"),
        ("ntr=0", "ntr=1"),
        ("nmropt=0", "nmropt=1"),
    ],
)
def test_input_drift_rejected_before_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, old: str, new: str
) -> None:
    campaign, start, engine, _, _ = setup_native(tmp_path, monkeypatch)
    with pytest.raises(ValueError):
        native.run_native(campaign, campaign / "bad", start, MDIN.replace(old, new))
    engine.assert_not_called()


def test_engine_and_scc_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign, start, _, _, settings = setup_native(tmp_path, monkeypatch)
    settings["exit"] = 1
    with pytest.raises(RuntimeError, match="engine"):
        native.run_native(campaign, campaign / "failed", start, MDIN)
    settings["exit"] = 0
    monkeypatch.setattr(native, "audit_dftb_scc", lambda text: {"passed": False})
    with pytest.raises(RuntimeError, match="SCC"):
        native.run_native(campaign, campaign / "scc", start, MDIN)
    assert (campaign / "scc/electronic_failure.json").exists()


def test_foreign_fixed_environment_rejected_before_engine(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject a starting restart with different frozen coordinates."""
    campaign, start, engine, _, _ = setup_native(tmp_path, monkeypatch)
    changed = np.zeros((2, 3))
    changed[1, 0] = 0.001
    monkeypatch.setattr(
        native.Rst7, "open", lambda path: SimpleNamespace(coordinates=changed)
    )
    with pytest.raises(ValueError, match="fixed environment"):
        native.run_native(campaign, campaign / "foreign", start, MDIN)
    engine.assert_not_called()
