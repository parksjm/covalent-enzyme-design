"""Optional integration checks against explicitly selected retained job inputs."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
from dataclasses import replace

import numpy as np
from parmed.amber import AmberParm
import pytest

from nyl12_covalent import local_repair_run
from nyl12_covalent import local_repair_model
from nyl12_covalent import local_water_premin
from nyl12_covalent.local_repair_checks import (
    audit_local_model,
    local_restraint_selection,
)
from nyl12_covalent.local_repair_model import identity, prepare_model

FIXTURE = os.environ.get("NYL12_REPAIR_FIXTURE")
pytestmark = pytest.mark.skipif(
    not FIXTURE, reason="Set NYL12_REPAIR_FIXTURE to retained preparation job directory"
)


@pytest.fixture
def paths() -> tuple[Path, Path, Path, Path]:
    root = Path(FIXTURE)
    return (
        root / "source.prmtop",
        root / "source.rst7",
        root / "af3.cif",
        root / "calculation",
    )


def test_exact_mapping_parameters_and_original_oxygen(
    paths: tuple, tmp_path: Path
) -> None:
    topology, coordinates, cif, _ = paths
    report = prepare_model(topology, coordinates, cif, tmp_path)
    assert report["water_seed"]["oxygen_displacement_A"] == 0
    assert report["qm_atoms"] == 109
    assert report["qm_electrons_with_link_H"] == 364
    assert report["water_parameters"][0]["epsilon_kcal_mol"] == pytest.approx(0.152)
    model = AmberParm(
        str(tmp_path / "model.prmtop"), xyz=str(tmp_path / "prepared.rst7")
    )
    audit = audit_local_model(model)
    assert audit["inventory_passed"] and audit["bond_screen_passed"]
    assert not audit["geometry_screen_passed"]
    selection = local_restraint_selection(model)
    assert {233, 234, 329, 330}.issubset(selection["local_residues_1based"])
    assert selection["all_hydrogens_mobile"]
    assert identity(232) == ("E", 1)
    with pytest.raises(ValueError, match="index"):
        identity(330)


@pytest.mark.parametrize("defect", ["proton", "water_name", "nonfinite", "ligand"])
def test_chemistry_defects_fail_closed(paths: tuple, defect: str) -> None:
    _, _, _, calculation = paths
    model = AmberParm(
        str(calculation / "model.prmtop"), xyz=str(calculation / "prepared.rst7")
    )
    if defect == "proton":
        next(a for a in model.residues[232].atoms if a.name == "H3").name = "HG1"
    elif defect == "water_name":
        model.residues[329].name = "HOH"
    elif defect == "ligand":
        model.residues[328].name = "BAD"
    else:
        coordinates = np.array(model.coordinates, copy=True)
        coordinates[0, 0] = np.nan
        model.coordinates = coordinates
    with pytest.raises(ValueError):
        audit_local_model(model)


def test_invalid_stage_input_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        local_repair_run.write_stage_input(tmp_path, ":1", [[1, 2]], 1, 0.05, 5)


def test_water_first_reorder_maps_back_and_preserves_fixed_atoms(
    paths: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, _, _, calculation = paths
    for name in ("model.prmtop", "prepared.rst7"):
        shutil.copy2(calculation / name, tmp_path / name)

    def fake_engine(command: list[str], **kwargs: object) -> dict:
        shutil.copy2(command[command.index("-c") + 1], command[command.index("-r") + 1])
        Path(command[command.index("-o") + 1]).write_text(
            "FINAL RESULTS\n NSTEP ENERGY RMS GMAX NAME NUMBER\n 1 -100.0 0.0001 0.001 O 1\n"
        )
        return {"exit_status": 0, "stopped_on_scf_failure": False}

    monkeypatch.setattr(local_water_premin, "run_checked_engine", fake_engine)
    restart = local_water_premin.preminimize_water(
        tmp_path, Path("/test/sander"), {}, [[3477, 5052], [3489, 5052]]
    )
    assert restart.name == "original_order.rst7"
    report = json.loads((tmp_path / "00_tip3p_water_only/acceptance.json").read_text())
    assert report["nonwater_maximum_displacement_A"] < 1e-6
    assert report["oxygen_displacement_from_AF3_A"] < 1e-6


@pytest.mark.parametrize("defect", ["duplicate", "residue", "extra", "missing"])
def test_invalid_af3_atom_mapping_fails(
    paths: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, defect: str
) -> None:
    topology, coordinates, cif, _ = paths
    sites = local_repair_model.read_atom_sites(cif)
    if defect == "duplicate":
        sites.append(sites[0])
    elif defect == "residue":
        sites[0] = replace(sites[0], residue="BAD")
    elif defect == "extra":
        sites.append(replace(sites[0], chain="Z"))
    else:
        sites.pop(0)
    monkeypatch.setattr(local_repair_model, "read_atom_sites", lambda path: sites)
    with pytest.raises(ValueError):
        prepare_model(topology, coordinates, cif, tmp_path)


@pytest.mark.parametrize("engine_exit", [0, 9])
def test_engine_exit_and_scientific_acceptance_are_separate(
    paths: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, engine_exit: int
) -> None:
    topology, coordinates, cif, _ = paths
    amber_bin = tmp_path / "amber/bin"
    amber_bin.mkdir(parents=True)
    (amber_bin / "sander").write_text("#!/bin/sh\nexit 0\n")
    (amber_bin / "sander").chmod(0o755)
    data = amber_bin.parent / "dat/slko"
    data.mkdir(parents=True)
    (data / "fixture.dat").write_text(
        "explicit test fixture, not a DFTB parameter set\n"
    )
    calls = []

    def fake_engine(command: list[str], **kwargs: object) -> dict:
        calls.append(command)
        if engine_exit:
            return {"exit_status": engine_exit, "stopped_on_scf_failure": False}
        shutil.copy2(command[command.index("-c") + 1], command[command.index("-r") + 1])
        Path(command[command.index("-o") + 1]).write_text(
            "FINAL RESULTS\n NSTEP ENERGY RMS GMAX NAME NUMBER\n 1 -100.0 0.01 0.02 O 1\n"
        )
        return {"exit_status": 0, "stopped_on_scf_failure": False}

    monkeypatch.setattr(local_repair_run, "run_checked_engine", fake_engine)
    monkeypatch.setattr(
        local_repair_run,
        "preminimize_water",
        lambda output, *args: output / "prepared.rst7",
    )
    destination = tmp_path / "result"
    if engine_exit:
        with pytest.raises(RuntimeError, match="Sander"):
            local_repair_run.run_repair(
                topology, coordinates, cif, amber_bin, destination
            )
        assert len(calls) == 1
        process = json.loads(
            (destination / "01_hydrogens_water/process.json").read_text()
        )
        assert process["exit_status"] == 9
        assert not (destination / "02_local_pocket").exists()
    else:
        local_repair_run.run_repair(topology, coordinates, cif, amber_bin, destination)
        assert len(calls) == 3
        result = json.loads((destination / "acceptance.json").read_text())
        assert not result["accepted_local_geometry"]
        assert not result["stationary_point_validated"]
        references = [c[c.index("-ref") + 1] for c in calls]
        assert len(set(references)) == 1
        released = (destination / "03_bridge_release/bridge.RST").read_text()
        assert "rk2=0.0" in released and "rk3=0.0" in released
