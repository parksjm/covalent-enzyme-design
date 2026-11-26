"""Focused scientific input and stage acceptance checks for hydration."""

from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from nyl12_covalent import deacylation_hydration as h

SOURCE = """&qmmm
qm_theory='DFTB3', qmcharge=-2, dftb_telec=0.0,
qmmask='@1-98', scfconv=1.d-12,
/
&dlfind
active=':1-2',
/
"""


def test_inputs() -> None:
    for stage in ("heat", "equil"):
        text = h.build_input(SOURCE, ":1-2", stage, 123)
        assert "qmmask='@1-98'" in text and "qmshake=0" in text
        assert "ibelly=1" in text and "ntr=0" in text and "ntc=1, ntf=1" in text
        assert "DISANG" not in text and "DUMPAVE" not in text
        assert ("TEMP0" in text) == (stage == "heat")
        assert ("nmropt=1" in text) == (stage == "heat")
    with pytest.raises(ValueError):
        h.build_input(SOURCE, ":2", "heat", 123)
    with pytest.raises(ValueError):
        h.build_input(SOURCE, ":1-2", "bad", 123)
    with pytest.raises(ValueError):
        h.build_input(SOURCE, ":1-2", "heat", 0)


def test_output_gates(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        h, "audit_dftb_scc", lambda text: {"verified_scc_history": True}
    )
    text = "NSTEP = 2000 TIME(PS) = 1.0 TEMP(K) = 300.0\n EPtot = -12.5\n"
    assert h.audit_output(text, 2000)["completed_steps"] == 2000
    with pytest.raises(ValueError, match="Temperature"):
        h.audit_output(text.replace("2000", "10000").replace("300.0", "200.0"), 10000)
    for broken in (
        text.replace("2000", "1999"),
        text.replace("-12.5", "NaN"),
        text.replace("300.0", "400.0"),
    ):
        with pytest.raises(ValueError):
            h.audit_output(broken, 2000)
    monkeypatch.setattr(
        h, "audit_dftb_scc", lambda text: {"verified_scc_history": False}
    )
    with pytest.raises(ValueError, match="SCC"):
        h.audit_output(text, 2000)


def test_restart_gates(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    base = np.zeros((2, 3))
    active = np.array([True, False])
    box = np.array([10, 10, 10, 90, 90, 90])
    state = SimpleNamespace(coordinates=base.copy(), vels=base.copy(), box=box.copy())
    monkeypatch.setattr(h.Rst7, "open", lambda path: state)
    _, result = h.audit_restart(tmp_path / "r", base, active, box)
    assert result["frozen_max_displacement_A"] == 0
    state.coordinates[1, 0] = 1e-3
    with pytest.raises(ValueError, match="exterior"):
        h.audit_restart(tmp_path / "r", base, active, box)
    state.coordinates = base.copy()
    state.box[0] = 11
    with pytest.raises(ValueError, match="box"):
        h.audit_restart(tmp_path / "r", base, active, box)
    state.box = box.copy()
    state.vels[0, 0] = np.nan
    with pytest.raises(ValueError, match="Nonfinite"):
        h.audit_restart(tmp_path / "r", base, active, box)


def test_prepared_runner_and_chemical_stop(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise preparation, frozen checks, stage chaining and chemistry stop."""
    import json
    from nyl12_covalent.deacylation_hydration_prepare import prepare

    repo = Path(__file__).resolve().parents[1]
    sources = {}
    for key in (
        "endpoint_restart",
        "source_topology",
        "endpoint_native_mdin",
        "sander",
    ):
        path = tmp_path / key
        path.write_text(SOURCE if key == "endpoint_native_mdin" else "input")
        sources[key] = {"path": str(path), "sha256": h.digest(path)}
    audit = tmp_path / "audit.json"
    payload = {
        "files": sources,
        "native_model": {
            "active_mask": ":1-2",
            "qmmask": "@1-98",
            "atom_count": 104922,
            "qm_atoms": 98,
            "qmcharge": -2,
            "link_atoms": 8,
            "active_atoms": 4165,
            "active_shell_A": 12,
        },
        "independent_acceptance": {
            "force_converged": True,
            "electronic_history_passed": True,
            "force_max_hartree_per_bohr": 9e-6,
        },
    }
    audit.write_text(json.dumps(payload))
    campaign = tmp_path / "campaign"
    prepare(campaign, audit, repo, 123)
    assert (
        campaign / "source_snapshot/nyl12_covalent/deacylation_hydration_prepare.py"
    ).exists()
    base = np.zeros((104922, 3))
    base[13915, 0] = 1.3
    active = np.zeros(104922, dtype=bool)
    active[:4165] = True
    box = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    qm = np.zeros(104922, dtype=bool)
    qm[:98] = True
    bonds = [
        SimpleNamespace(
            atom1=SimpleNamespace(idx=i), atom2=SimpleNamespace(idx=100 + i)
        )
        for i in range(8)
    ]
    monkeypatch.setattr(h, "AmberParm", lambda p: SimpleNamespace(bonds=bonds))
    monkeypatch.setattr(
        h,
        "AmberMask",
        lambda parm, mask: SimpleNamespace(
            Selection=lambda: qm if mask == "@1-98" else active
        ),
    )
    monkeypatch.setattr(
        h.Rst7, "open", lambda p: SimpleNamespace(coordinates=base, box=box)
    )
    monkeypatch.setattr(h, "_journal", lambda *args: None)
    monkeypatch.setattr(
        h, "audit_output", lambda text, steps: {"completed_steps": steps}
    )
    monkeypatch.setattr(
        h,
        "audit_restart",
        lambda path, baseline, mask, box: (base, {"box_unchanged": True}),
    )

    def fake_run(cmd: list[str], **kwargs: object) -> SimpleNamespace:
        Path(cmd[cmd.index("-o") + 1]).write_text("mdout")
        Path(cmd[cmd.index("-r") + 1]).write_text("restart")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(h.subprocess, "run", fake_run)
    assert h.run(campaign)["status"] == "completed"
    assert json.loads((campaign / "equil/scientific_acceptance.json").read_text())[
        "chemistry"
    ]["acyl_bond_retained"]
    campaign2 = tmp_path / "campaign2"
    prepare(campaign2, audit, repo, 123)
    base[13915, 0] = 2.1
    assert h.run(campaign2)["reason"] == "acyl_bond_not_retained"
    assert not (campaign2 / "equil").exists()
    sources["source_topology"]["sha256"] = "wrong"
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="hash"):
        prepare(tmp_path / "bad", audit, repo, 123)
    payload["independent_acceptance"]["force_converged"] = False
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="force/SCC"):
        prepare(tmp_path / "unaccepted", audit, repo, 123)
