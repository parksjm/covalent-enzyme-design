"""Provenance and scientific acceptance checks for Gaussian conformer pairs."""

from __future__ import annotations

import json
import shlex
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from nyl12_covalent.gaussian_checks import assess_gaussian, final_orientation
from nyl12_covalent.gaussian_conformer_geometry import (
    dihedral_degrees,
    named_psi_indices,
    parse_stage_energies,
    sha256,
)

REQUIRED_ACCEPTANCE = (
    "accepted",
    "completed",
    "optimization_completed",
    "frequency_minimum_passed",
    "electronic_state_passed",
    "esp_grid_present",
)


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object and reject non-object data."""
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _verify_hash(path: Path, expected: str, label: str) -> str:
    """Verify one recorded file hash and return the observed digest."""
    if not path.is_file():
        raise ValueError(f"Recorded {label} file is missing: {path}")
    observed = sha256(path)
    if observed != expected:
        raise ValueError(f"{label} SHA-256 mismatch: {path}")
    return observed


def verify_provenance(
    conformer: str, prepared_dir: Path, calculation_dir: Path
) -> dict[str, Any]:
    """Verify input, source, executable, and process records before analysis."""
    calculation = read_json(calculation_dir / "calculation.json")
    identity = calculation.get("identity", {})
    expected_inputs = identity.get("input_sha256", {})
    verified_inputs: dict[str, str] = {}
    for name in ("input.gjf", "molecule.sdf", "atom_map.json"):
        if name not in expected_inputs:
            raise ValueError(f"Missing recorded input hash for {name}: {conformer}")
        for folder in (prepared_dir, calculation_dir):
            verified_inputs[f"{folder.name}/{name}"] = _verify_hash(
                folder / name, expected_inputs[name], f"input {name}"
            )
    if calculation.get("source") != str(prepared_dir.resolve()):
        raise ValueError(f"Prepared source path mismatch: {conformer}")

    source_hashes = calculation.get("source_sha256", {})
    if not source_hashes:
        raise ValueError(f"No source snapshot hashes recorded: {conformer}")
    verified_source = {}
    for name, expected in source_hashes.items():
        verified_source[name] = _verify_hash(
            calculation_dir / "source" / name, expected, f"source {name}"
        )
    gaussian_source = source_hashes.get("gaussian_checks.py")
    if gaussian_source != identity.get("acceptance_source_sha256"):
        raise ValueError(f"Acceptance source hash records disagree: {conformer}")
    project_root = Path(__file__).resolve().parents[2]
    if (
        sha256(project_root / "src/nyl12_covalent/gaussian_checks.py")
        != gaussian_source
    ):
        raise ValueError(
            "Current Gaussian acceptance code differs from calculation snapshot"
        )
    shell = calculation_dir / "run_gaussian.sh"
    _verify_hash(shell, calculation.get("shell_sha256", ""), "Gaussian shell")

    executable_lines = [
        line.strip()
        for line in shell.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    executable = next(
        (Path(shlex.split(line)[0]) for line in executable_lines if "/g16/g16" in line),
        None,
    )
    profile_line = next(
        (line for line in executable_lines if line.startswith("source ")), None
    )
    esp_command = read_json(calculation_dir / "espgen_command.json").get("argv", [])
    if executable is None or profile_line is None or not esp_command:
        raise ValueError(f"Cannot identify Gaussian/ESP executables: {conformer}")
    engine_paths = {
        "g16": executable,
        "profile": Path(shlex.split(profile_line)[1]),
        "espgen": Path(esp_command[0]),
    }
    verified_engines = {}
    for name, path in engine_paths.items():
        expected = identity.get("engine_sha256", {}).get(name)
        if not expected:
            raise ValueError(f"Missing recorded {name} hash: {conformer}")
        verified_engines[name] = _verify_hash(path, expected, f"engine {name}")

    process = read_json(calculation_dir / "process.json")
    if (
        process.get("exit_status") != 0
        or not process.get("started_utc")
        or not process.get("ended_utc")
    ):
        raise ValueError(f"Gaussian process failed or is incomplete: {conformer}")
    return {
        "input_sha256": verified_inputs,
        "source_sha256": verified_source,
        "engine_sha256": verified_engines,
        "shell_sha256": calculation["shell_sha256"],
        "process_exit_status": process["exit_status"],
    }


def validate_acceptance_record(acceptance: dict[str, Any], conformer: str) -> None:
    """Require every recorded QM prerequisite and stereochemical gate to pass."""
    for field in REQUIRED_ACCEPTANCE:
        if acceptance.get(field) is not True:
            raise ValueError(
                f"Required acceptance field {field} is not true: {conformer}"
            )
    if acceptance.get("fitted_forcefield_ready") is not False:
        raise ValueError(f"Unexpected force-field readiness claim: {conformer}")
    geometry = acceptance.get("geometry", {})
    if geometry.get("accepted") is not True or geometry.get(
        "expected_stereocenters"
    ) != geometry.get("observed_stereocenters"):
        raise ValueError(f"Recorded symmetry/chirality gate failed: {conformer}")


def audit_conformer(
    conformer: str, prepared_dir: Path, calculation_dir: Path
) -> tuple[dict[str, Any], np.ndarray, list[str]]:
    """Verify records and independently reassess one conformer's Gaussian output."""
    provenance = verify_provenance(conformer, prepared_dir, calculation_dir)
    acceptance = read_json(calculation_dir / "acceptance.json")
    validate_acceptance_record(acceptance, conformer)
    geometry = acceptance["geometry"]
    log_path = calculation_dir / "gaussian.log"
    molecule_path = calculation_dir / "molecule.sdf"
    independently_assessed = assess_gaussian(log_path, molecule_path)
    if not independently_assessed.get("accepted"):
        raise ValueError(f"Independent QM assessment failed: {conformer}")
    recorded_frequencies = acceptance.get("frequencies_cm_inverse", [])
    assessed_frequencies = independently_assessed["frequencies_cm_inverse"]
    if (
        acceptance.get("expected_frequency_count") != 57
        or len(recorded_frequencies) != 57
        or not np.allclose(
            recorded_frequencies, assessed_frequencies, rtol=0, atol=1e-8
        )
    ):
        raise ValueError(
            f"The accepted record does not match all 57 Gaussian modes: {conformer}"
        )
    text = log_path.read_text(encoding="utf-8", errors="replace")
    numbers, xyz = final_orientation(text)
    recorded_xyz = np.asarray(acceptance.get("final_coordinates_A"), dtype=float)
    if recorded_xyz.shape != xyz.shape or not np.allclose(
        recorded_xyz, xyz, atol=1e-6, rtol=0
    ):
        raise ValueError(
            f"Accepted coordinates disagree with Gaussian log: {conformer}"
        )
    atoms = read_json(calculation_dir / "atom_map.json").get("atoms", [])
    if len(atoms) != len(numbers) or [atom["index"] for atom in atoms] != list(
        range(len(atoms))
    ):
        raise ValueError(f"Atom map does not match Gaussian atom order: {conformer}")
    elements = [atom["element"] for atom in atoms]
    if [
        Chem.GetPeriodicTable().GetAtomicNumber(element) for element in elements
    ] != numbers:
        raise ValueError(f"Atom map elements disagree with Gaussian order: {conformer}")
    energy = parse_stage_energies(text)
    result = {
        "name": conformer,
        "optimized_wb97xd_energy_Eh": energy["optimized_wb97xd_Eh"],
        "hf_esp_energy_Eh": energy["hf_esp_Eh"],
        "psi_N_CA_C_NME_N_degrees": dihedral_degrees(xyz[named_psi_indices(atoms)]),
        "expected_stereocenters": geometry["expected_stereocenters"],
        "observed_stereocenters": geometry["observed_stereocenters"],
        "frequency_count": len(acceptance["frequencies_cm_inverse"]),
        "frequencies_cm_inverse": recorded_frequencies,
        "minimum_frequency_cm_inverse": min(recorded_frequencies),
        "artifact_sha256": {
            name: sha256(calculation_dir / name)
            for name in ("gaussian.log", "acceptance.json", "process.json", "esp.dat")
        },
        "minimum_passed": independently_assessed["frequency_minimum_passed"],
        "scientific_acceptance_recomputed": independently_assessed["accepted"],
        "provenance": provenance,
        "optimized_coordinates_A": xyz.tolist(),
        "atom_map": atoms,
        "prepared_molecule_sha256": sha256(prepared_dir / "molecule.sdf"),
    }
    return result, xyz, elements
