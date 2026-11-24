"""Tests for scientific progress publication and source validation."""

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from nyl12_covalent import nterminal_relay_progress as progress


def fixture_files(tmp_path: Path) -> tuple[Path, Path, Path, list[dict]]:
    """Prepare a minimal mocked source and two progress records."""
    campaign = tmp_path / "campaign"
    model = campaign / "model/inputs"
    template = campaign / "template"
    model.mkdir(parents=True)
    template.mkdir()
    topology = model / "system.prmtop"
    topology.write_text("topology")
    (template / "stripped.parm7").write_text("stripped")
    (template / "view_relay.vmd").write_text("mol addfile optimization.dcd\n")
    (template / "manifest.json").write_text(
        json.dumps(
            {
                "topology_sha256": progress.sha(topology),
                "stripped_to_original_atom_indices_0based": [0, 2],
            }
        )
    )
    cpptraj = tmp_path / "cpptraj"
    cpptraj.write_text("program")
    (campaign / "config.json").write_text(
        json.dumps(
            {
                "model": "model",
                "viewer_template": "template",
                "cpptraj": str(cpptraj),
                "source_topology_sha256": progress.sha(topology),
                "reference_energy_hartree": -100.0,
            }
        )
    )
    distances = {
        f"{h}_{a}": 1.0 for h in ("H1", "H2", "HG", "HN") for a in ("N", "OG", "NL")
    }
    distances.update(C_OG=1.5, C_NL=1.5, C_O=1.2)
    rows = []
    for i in range(2):
        restart = tmp_path / f"restart{i}"
        restart.write_text("restart")
        d = distances.copy()
        d["HG_NL"] += i * 0.4
        rows.append(
            {
                "label": "window033_start" if i == 0 else "window_001",
                "restart": str(restart),
                "accepted": True,
                "native_converged": True,
                "geometry": {"distances_A": d},
                "total_kcal_mol": -10 + i,
                "bias_kcal_mol": 0.0,
                "physical_kcal_mol": -10 + i,
                "target_A": float(i),
                "observed_coordinate_A": float(i),
            }
        )
    return campaign, template, cpptraj, rows


def setup_mocks(monkeypatch: pytest.MonkeyPatch, template: Path) -> None:
    """Mock only native readers and external conversion, retaining publisher logic."""
    atoms = [
        SimpleNamespace(idx=i, residue=SimpleNamespace(name=n))
        for i, n in enumerate(("ALA", "WAT", "P6I"))
    ]
    monkeypatch.setattr(
        progress.parmed,
        "load_file",
        lambda path: SimpleNamespace(
            atoms=atoms
            if not str(path).endswith("stripped.parm7")
            else [atoms[0], atoms[2]]
        ),
    )
    monkeypatch.setattr(
        progress.parmed.amber.Rst7,
        "open",
        lambda path: SimpleNamespace(
            coordinates=np.arange(9).reshape(3, 3),
            box=np.array([10, 10, 10, 90, 90, 90]),
        ),
    )

    def run(args: list[str], **kwargs: object) -> SimpleNamespace:
        script = Path(args[-1]).read_text()
        output = script.split('trajout "')[1].split('"')[0]
        Path(output).write_bytes(b"mock-dcd")
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(progress.subprocess, "run", run)
    monkeypatch.setattr(
        progress,
        "validate_dcd",
        lambda path, states: {"frames": len(states), "atoms": len(states[0][0])},
    )


def test_publish(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign, template, _, rows = fixture_files(tmp_path)
    setup_mocks(monkeypatch, template)
    out = tmp_path / "branch"
    progress.publish(campaign, out, rows)
    doc = json.loads((out / "progress.json").read_text())
    assert doc["dcd_validation"] == {"frames": 2, "atoms": 2}
    assert doc["distance_changes"][1]["large_signed_changes_A"][
        "HG_NL"
    ] == pytest.approx(0.4)
    assert "reaction_path.dcd" in (out / "vmd_progress/view_path.vmd").read_text()
    assert (
        "strip :WAT,Na+,Cl-,NA,CL" in (out / "vmd_progress/export.cpptraj").read_text()
    )
    assert (out / "plots/energies.png").exists() and (
        out / "plots/distances.png"
    ).exists()
    assert (out / "progress.csv").read_text().count("\n") == 3


def test_failures(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign, template, cpptraj, rows = fixture_files(tmp_path)
    setup_mocks(monkeypatch, template)
    with pytest.raises(ValueError):
        progress.publish(campaign, tmp_path / "out", [])
    cfg = json.loads((campaign / "config.json").read_text())
    cfg["source_topology_sha256"] = "wrong"
    (campaign / "config.json").write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="hash"):
        progress.publish(campaign, tmp_path / "out", rows)
    cfg["source_topology_sha256"] = progress.sha(
        campaign / "model/inputs/system.prmtop"
    )
    (campaign / "config.json").write_text(json.dumps(cfg))
    monkeypatch.setattr(
        progress.subprocess,
        "run",
        Mock(return_value=SimpleNamespace(returncode=1, stdout="", stderr="failure")),
    )
    with pytest.raises(RuntimeError, match="cpptraj"):
        progress.publish(campaign, tmp_path / "out", rows)
