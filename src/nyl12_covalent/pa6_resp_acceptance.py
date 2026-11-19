"""Gate optimized PA6 geometries before fixed-geometry HF/MK ESP jobs."""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
import re
from typing import Any

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.gaussian_checks import final_orientation, validate_geometry
from nyl12_covalent.gaussian_conformer_geometry import dihedral_degrees, sha256
from nyl12_covalent.pa6_resp_esp import hf_esp_input
from nyl12_covalent.ti1_reference_graph import validate_pa6_zwitterion

LOGGER = logging.getLogger(__name__)
N_H_BOND_RANGE_ANGSTROM = (0.80, 1.25)
CARBOXYLATE_O_H_MIN_ANGSTROM = 1.20


def prepare_pa6_hf_esp(
    source_mol2: Path,
    preparation_manifest: Path,
    conformer_label: str,
    gaussian_log: Path,
    torsion_tolerance_degrees: float,
    acceptance_path: Path,
    esp_input_path: Path,
) -> dict[str, Any]:
    """Assess one DFT optimization and write an HF/MK input only if accepted.

    Args:
        source_mol2: Canonical 79-atom PA6 L4 zwitterion MOL2.
        preparation_manifest: Manifest from the constrained-input preparation.
        conformer_label: Exact conformer label in that manifest.
        gaussian_log: Completed DFT geometry-optimization output.
        torsion_tolerance_degrees: Maximum circular drift for each frozen torsion.
        acceptance_path: New ``scientific_acceptance.json`` path.
        esp_input_path: New fixed-geometry HF/MK Gaussian input path.

    Returns:
        Scientific acceptance record. A rejected result is written before a
        ``RuntimeError`` is raised.

    Raises:
        FileExistsError: If either output already exists.
        ValueError: If thresholds or declared provenance inputs are invalid.
        RuntimeError: If the optimized geometry fails scientific acceptance.
    """
    if acceptance_path.exists() or esp_input_path.exists():
        raise FileExistsError("Acceptance and ESP outputs must both be new")
    if not math.isfinite(torsion_tolerance_degrees) or not (
        0 < torsion_tolerance_degrees <= 30
    ):
        raise ValueError("Torsion tolerance must be finite and in (0, 30] degrees")
    try:
        report, elements, xyz = _assess_geometry(
            source_mol2,
            preparation_manifest,
            conformer_label,
            gaussian_log,
            torsion_tolerance_degrees,
        )
    except (OSError, ValueError, json.JSONDecodeError) as error:
        report = {
            "schema_version": 1,
            "accepted": False,
            "conformer": conformer_label,
            "failure": str(error),
            "esp_input_written": False,
        }
        acceptance_path.parent.mkdir(parents=True, exist_ok=True)
        acceptance_path.write_text(json.dumps(report, indent=2) + "\n")
        raise RuntimeError(f"PA6 geometry assessment failed: {error}") from error
    acceptance_path.parent.mkdir(parents=True, exist_ok=True)
    acceptance_path.write_text(json.dumps(report, indent=2) + "\n")
    if not report["accepted"]:
        raise RuntimeError("PA6 geometry failed scientific acceptance")
    esp_input_path.parent.mkdir(parents=True, exist_ok=True)
    esp_input_path.write_text(hf_esp_input(elements, xyz, conformer_label))
    report["esp_input_written"] = True
    report["esp_input"] = str(esp_input_path.resolve())
    report["esp_input_sha256"] = sha256(esp_input_path)
    acceptance_path.write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info("Accepted %s and prepared fixed-geometry HF/MK input", conformer_label)
    return report


def _assess_geometry(
    source_mol2: Path,
    manifest_path: Path,
    label: str,
    log_path: Path,
    tolerance: float,
) -> tuple[dict[str, Any], list[int], NDArray[np.float64]]:
    substrate = load_pa6(source_mol2)
    validate_pa6_zwitterion(substrate.mol)
    names = [atom.GetProp("name") for atom in substrate.mol.GetAtoms()]
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("source_mol2_sha256") != sha256(source_mol2):
        raise ValueError("Source MOL2 hash differs from preparation manifest")
    if manifest.get("atom_names_in_canonical_order") != names:
        raise ValueError("Manifest atom names/order differ from canonical MOL2")
    records = [
        item for item in manifest.get("conformers", []) if item.get("label") == label
    ]
    if len(records) != 1:
        raise ValueError(f"Expected one manifest conformer labeled {label!r}")
    record = records[0]
    initial_path = manifest_path.parent / record["coordinates_xyz"]
    if sha256(initial_path) != record.get("coordinates_xyz_sha256"):
        raise ValueError("Prepared conformer XYZ hash mismatch")
    initial_elements, initial_xyz = _read_xyz(initial_path)
    canonical_elements = [atom.GetAtomicNum() for atom in substrate.mol.GetAtoms()]
    if initial_elements != canonical_elements:
        raise ValueError("Prepared XYZ elements/order differ from canonical MOL2")

    text = log_path.read_text(errors="replace")
    elements, xyz = final_orientation(text)
    geometry = validate_geometry(substrate.mol, elements, xyz)
    normal_termination = (
        text.count("Normal termination of Gaussian 16") == 1
        and "Normal termination of Gaussian 16" in text.rstrip().splitlines()[-1]
        and "Error termination" not in text
    )
    charges = [
        (int(charge), int(spin))
        for charge, spin in re.findall(
            r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", text
        )
    ]
    expected_method, expected_basis, expected_solvent = _expected_method_metadata(
        manifest
    )
    route_lines = [
        line.lower()
        for line in text.splitlines()
        if "#p" in line.lower() and re.search(r"\w+/\w+", line)
    ]
    route_methods = {
        tuple(part.lower() for part in match.groups())
        for line in route_lines
        if (match := re.search(r"([a-z0-9-]+)/([a-z0-9()+_-]+)", line))
    }
    route_solvents = {
        match.group(1).lower()
        for line in route_lines
        if (match := re.search(r"scrf=\(smd,solvent=([a-z0-9_-]+)\)", line))
    }
    if len(route_methods) > 1 or len(route_solvents) > 1:
        raise ValueError("Gaussian log contains ambiguous method or solvent routes")
    method = route_methods == {
        (expected_method.lower(), expected_basis.lower())
    } and route_solvents == {expected_solvent.lower()}
    electronic_state = bool(charges) and all(item == (0, 1) for item in charges)
    torsions = _torsion_report(
        manifest.get("frozen_torsions"), names, initial_xyz, xyz, tolerance
    )
    protonation = _terminal_protonation(substrate.mol, names, xyz)
    checks = {
        "normal_termination": normal_termination,
        "optimization_completed": "Optimization completed." in text,
        "method_basis_and_smd_solvent_match_manifest": method,
        "electronic_state_0_1": electronic_state,
        "canonical_elements_and_connectivity": bool(geometry["accepted"]),
        "frozen_torsions_within_tolerance": all(item["passed"] for item in torsions),
        "terminal_proton_inventory_retained": protonation["passed"],
    }
    accepted = all(checks.values())
    report = {
        "schema_version": 1,
        "accepted": accepted,
        "scope": "PA6 L4 constrained DFT geometry acceptance before RESP ESP",
        "conformer": label,
        "checks": checks,
        "thresholds": {
            "torsion_tolerance_degrees": tolerance,
            "n_terminal_n_h_bond_range_angstrom": list(N_H_BOND_RANGE_ANGSTROM),
            "carboxylate_o_h_minimum_angstrom": CARBOXYLATE_O_H_MIN_ANGSTROM,
        },
        "torsions": torsions,
        "terminal_protonation": protonation,
        "geometry": geometry,
        "source_sha256": {
            "mol2": sha256(source_mol2),
            "preparation_manifest": sha256(manifest_path),
            "initial_xyz": sha256(initial_path),
            "gaussian_log": sha256(log_path),
        },
        "final_coordinates_angstrom": xyz.tolist(),
        "esp_input_written": False,
    }
    return report, elements, xyz


def _expected_method_metadata(manifest: dict[str, Any]) -> tuple[str, str, str]:
    """Read explicit restart settings or unambiguously parse a legacy route."""
    gaussian = manifest.get("gaussian")
    if not isinstance(gaussian, dict):
        raise ValueError("Manifest lacks Gaussian method metadata")
    method, basis = gaussian.get("method"), gaussian.get("basis")
    solvent = gaussian.get("solvent")
    model = gaussian.get("solvation_model")
    if (
        method is not None
        or basis is not None
        or solvent is not None
        or model is not None
    ):
        if not all(
            isinstance(value, str) and value for value in (method, basis, solvent)
        ):
            raise ValueError(
                "Manifest method, basis, and solvent metadata are required"
            )
        if model != "SMD":
            raise ValueError("Manifest solvent model must explicitly be SMD")
        return method, basis, solvent
    route = gaussian.get("route")
    if not isinstance(route, str):
        raise ValueError("Manifest lacks explicit method/basis/solvent metadata")
    methods = re.findall(r"([A-Za-z0-9-]+)/([A-Za-z0-9()+_-]+)", route)
    solvents = re.findall(r"SCRF=\(SMD,Solvent=([A-Za-z0-9_-]+)\)", route, re.I)
    if len(methods) != 1 or len(solvents) != 1:
        raise ValueError(
            "Legacy Gaussian route lacks unambiguous method/basis/SMD solvent"
        )
    return methods[0][0], methods[0][1], solvents[0]


def _read_xyz(path: Path) -> tuple[list[int], NDArray[np.float64]]:
    lines = path.read_text().splitlines()
    if len(lines) < 2 or int(lines[0]) != len(lines) - 2:
        raise ValueError(f"Malformed XYZ: {path}")
    periodic = Chem.GetPeriodicTable()
    elements: list[int] = []
    coordinates = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"Malformed XYZ coordinate: {line}")
        elements.append(periodic.GetAtomicNumber(fields[0]))
        coordinates.append([float(value) for value in fields[1:]])
    xyz = np.asarray(coordinates, dtype=float)
    if not np.isfinite(xyz).all():
        raise ValueError("XYZ contains nonfinite coordinates")
    return elements, xyz


def _torsion_report(
    specification: object,
    names: list[str],
    initial: NDArray[np.float64],
    final: NDArray[np.float64],
    tolerance: float,
) -> list[dict[str, Any]]:
    if not isinstance(specification, list) or not specification:
        raise ValueError("Manifest lacks frozen torsion specifications")
    indices = {name: index for index, name in enumerate(names)}
    results = []
    for item in specification:
        atoms = item.get("atoms") if isinstance(item, dict) else None
        if (
            not isinstance(atoms, list)
            or len(atoms) != 4
            or any(atom not in indices for atom in atoms)
        ):
            raise ValueError("Invalid frozen torsion record in manifest")
        atom_indices = [indices[atom] for atom in atoms]
        start = dihedral_degrees(initial[atom_indices])
        end = dihedral_degrees(final[atom_indices])
        drift = abs((end - start + 180.0) % 360.0 - 180.0)
        results.append(
            {
                "name": item.get("name"),
                "atoms": atoms,
                "initial_degrees": start,
                "final_degrees": end,
                "absolute_circular_drift_degrees": drift,
                "passed": drift <= tolerance,
            }
        )
    return results


def _terminal_protonation(
    molecule: Chem.Mol, names: list[str], xyz: NDArray[np.float64]
) -> dict[str, Any]:
    by_name = {name: index for index, name in enumerate(names)}
    n_index = by_name["N56"]
    terminal_hydrogens = sorted(
        names[neighbor.GetIdx()]
        for neighbor in molecule.GetAtomWithIdx(n_index).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    )
    if terminal_hydrogens != ["H1", "H2", "H3"]:
        raise ValueError("Canonical N56 terminal-hydrogen inventory changed")
    heavy_indices = [
        atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1
    ]
    n_h = {}
    owners = {}
    for hydrogen_name in terminal_hydrogens:
        h_index = by_name[hydrogen_name]
        n_h[hydrogen_name] = float(np.linalg.norm(xyz[n_index] - xyz[h_index]))
        nearest = min(
            heavy_indices,
            key=lambda index: float(np.linalg.norm(xyz[index] - xyz[h_index])),
        )
        owners[hydrogen_name] = names[nearest]
    carboxylate = {
        oxygen: min(
            float(np.linalg.norm(xyz[by_name[oxygen]] - xyz[atom.GetIdx()]))
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() == 1
        )
        for oxygen in ("O48", "O51")
    }
    passed = (
        all(
            N_H_BOND_RANGE_ANGSTROM[0] <= value <= N_H_BOND_RANGE_ANGSTROM[1]
            for value in n_h.values()
        )
        and all(owner == "N56" for owner in owners.values())
        and all(value >= CARBOXYLATE_O_H_MIN_ANGSTROM for value in carboxylate.values())
    )
    return {
        "passed": passed,
        "n56_hydrogens": terminal_hydrogens,
        "n56_h_distances_angstrom": n_h,
        "terminal_h_nearest_heavy_atom": owners,
        "minimum_carboxylate_oxygen_hydrogen_distances_angstrom": carboxylate,
    }
