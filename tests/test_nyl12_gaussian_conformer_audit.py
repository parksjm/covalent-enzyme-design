"""Focused checks for staged-energy parsing and optimized-coordinate exports."""

import numpy as np
import pytest
from rdkit import Chem

from nyl12_covalent.gaussian_conformer_audit import validate_acceptance_record
from nyl12_covalent.gaussian_conformer_geometry import (
    heavy_atom_kabsch_rmsd,
    parse_stage_energies,
    write_optimized_sdf,
)


def test_stage_energy_parser_uses_dft_opt_and_separate_hf_esp() -> None:
    log = """SCF Done: E(RwB97XD) = -10.0 A.U.
SCF Done: E(RwB97XD) = -11.0 A.U.
Optimization completed.
Normal termination of Gaussian 16
SCF Done: E(RwB97XD) = -11.0 A.U.
Frequency calculation
Normal termination of Gaussian 16
SCF Done: E(RHF) = -8.0 A.U.
Normal termination of Gaussian 16
"""
    assert parse_stage_energies(log) == {
        "optimized_wb97xd_Eh": -11.0,
        "hf_esp_Eh": -8.0,
    }


@pytest.mark.parametrize(
    "log",
    [
        "SCF Done: E(RwB97XD) = -1.0 A.U.",
        "Normal termination of Gaussian 16\nNormal termination of Gaussian 16",
        """SCF Done: E(RwB97XD) = -1.0 A.U.
Normal termination of Gaussian 16
Normal termination of Gaussian 16
""",
    ],
)
def test_stage_energy_parser_rejects_missing_stage_energies(log: str) -> None:
    with pytest.raises(ValueError):
        parse_stage_energies(log)


def test_heavy_atom_kabsch_rmsd_removes_rigid_motion() -> None:
    first = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [8.0, 8.0, 8.0]]
    )
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    second = first @ rotation + np.array([4.0, -3.0, 2.0])
    assert heavy_atom_kabsch_rmsd(first, second, ["C", "N", "O", "H"]) < 1e-12


def test_sdf_export_uses_final_optimized_coordinates(tmp_path) -> None:
    molecule = Chem.AddHs(Chem.MolFromSmiles("O"))
    molecule_path = tmp_path / "input.sdf"
    with Chem.SDWriter(str(molecule_path)) as writer:
        writer.write(molecule)
    optimized = np.array([[1.0, 2.0, 3.0], [1.9, 2.0, 3.0], [0.8, 2.8, 3.0]])
    output = tmp_path / "optimized.sdf"

    write_optimized_sdf(molecule_path, optimized, output)

    exported = next(iter(Chem.SDMolSupplier(str(output), removeHs=False)))
    assert exported is not None
    assert exported.GetNumAtoms() == molecule.GetNumAtoms()
    assert exported.GetNumBonds() == molecule.GetNumBonds()
    assert np.allclose(exported.GetConformer().GetPositions(), optimized)
    assert (
        exported.GetProp("coordinate_source") == "final Gaussian optimized orientation"
    )


def test_acceptance_requires_all_gates_and_preserved_chirality() -> None:
    record = {
        "accepted": True,
        "completed": True,
        "optimization_completed": True,
        "frequency_minimum_passed": True,
        "electronic_state_passed": True,
        "esp_grid_present": True,
        "fitted_forcefield_ready": False,
        "geometry": {
            "accepted": True,
            "expected_stereocenters": {"1": "S", "2": "R"},
            "observed_stereocenters": {"1": "S", "2": "R"},
        },
    }
    validate_acceptance_record(record, "conformer_01")
    record["geometry"]["observed_stereocenters"] = {"1": "R", "2": "S"}
    with pytest.raises(ValueError, match="chirality"):
        validate_acceptance_record(record, "conformer_01")


def test_last_completed_optimization_is_used() -> None:
    log = """SCF Done: E(RwB97XD) = -10.0 A.U.
Optimization completed.
Normal termination of Gaussian 16
SCF Done: E(RwB97XD) = -11.5 A.U.
Optimization completed.
Normal termination of Gaussian 16
SCF Done: E(RHF) = -8.0 A.U.
Normal termination of Gaussian 16
"""
    assert parse_stage_energies(log)["optimized_wb97xd_Eh"] == -11.5
