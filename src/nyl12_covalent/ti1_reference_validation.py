"""Validate constrained TI1 QM geometry and Hessian prerequisites."""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdchem

from .gaussian_checks import final_orientation, validate_geometry
from .gaussian_fchk import parse_gaussian_fchk
from .ti1_reference_coordinates import TORSION_LABELS

LOGGER = logging.getLogger(__name__)
_EXPECTED_REFERENCE = {"chain": "A", "residue_number": 233}
_EXPECTED_PSI = 134.110276
_EXPECTED_CHI1 = -58.964418
_CAP_PHI_CHOICES = {-120.0, -60.0}
_TORSION_TOLERANCE = 0.10
_BOND_WINDOWS = {
    "THR:OG1": (1.30, 1.75),
    "TIL:O": (1.25, 1.65),
    "TIL:N": (1.25, 1.70),
    "TIL:acylCH3": (1.35, 1.75),
}
_REQUIRED_ORIGINS = {
    "THR:CA",
    "THR:CB",
    "THR:N",
    "THR:OG1",
    "THR:C",
    "ACE:C",
    "NME:N",
    "TIL:C*",
    "TIL:O",
    "TIL:N",
    "TIL:acylCH3",
}


def assess_reference_optimization(directory: Path) -> dict[str, Any]:
    """Assess a constrained DFT geometry as a scoped TI1 reference.

    Args:
        directory: Directory containing atom_map.json, molecule.sdf and
            optimization.log from one prepared capped-surrogate job.

    Returns:
        Geometry, torsion, electronic-state and reactive-bond checks. An
        accepted report is only a constrained QM geometry prerequisite.

    Raises:
        ValueError: If input metadata, molecular identity, or Gaussian
            orientation is missing or malformed.
        FileNotFoundError: If any required input file is absent.
    """
    metadata = _read_json(directory / "atom_map.json")
    molecule = _read_one_sdf(directory / "molecule.sdf")
    atom_map = _validate_model_identity(molecule, metadata)
    _validate_reference_contract(metadata, atom_map)

    text = (directory / "optimization.log").read_text(encoding="utf-8")
    atomic_numbers, coordinates = final_orientation(text)
    geometry = validate_geometry(molecule, atomic_numbers, coordinates)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, point in enumerate(coordinates):
        conformer.SetAtomPosition(index, tuple(float(value) for value in point))

    torsions = _measure_constraints(conformer, metadata["constraints"])
    state_records = [
        (int(charge), int(multiplicity))
        for charge, multiplicity in re.findall(
            r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", text
        )
    ]
    state_pass = bool(state_records) and all(
        charge == -1 and multiplicity == 1 for charge, multiplicity in state_records
    )
    complete = (
        _ends_with_normal_termination(text)
        and "Optimization completed." in text
        and "Error termination" not in text
    )
    distances = _reactive_bond_distances(atom_map, coordinates)
    bond_pass = all(
        lower <= distances[name] <= upper
        for name, (lower, upper) in _BOND_WINDOWS.items()
    )
    close_hydrogens = _oxyanion_hydrogen_contacts(
        metadata["atoms"], atom_map, coordinates
    )
    proton_transfer = any(contact["distance_A"] < 1.20 for contact in close_hydrogens)
    constraints_pass = all(
        row["error_degrees"] <= _TORSION_TOLERANCE for row in torsions
    )
    accepted = bool(
        complete
        and state_pass
        and geometry["accepted"]
        and constraints_pass
        and bond_pass
        and not proton_transfer
    )
    report = {
        "accepted": accepted,
        "scope": "constrained capped QM geometry prerequisite only",
        "forcefield_validated": False,
        "native_charge_transfer_eligible": False,
        "unconstrained_minimum_claimed": False,
        "optimization_completed_normally": complete,
        "electronic_state_passed": state_pass,
        "geometry": geometry,
        "constraints": torsions,
        "constraint_tolerance_degrees": _TORSION_TOLERANCE,
        "reactive_bonds_A": distances,
        "reactive_bond_windows_A": _BOND_WINDOWS,
        "reactive_bond_screen_passed": bond_pass,
        "oxyanion_close_hydrogens": close_hydrogens,
        "oxyanion_proton_transfer_detected": proton_transfer,
        "final_coordinates_A": coordinates.tolist(),
    }
    LOGGER.info("TI1 constrained reference geometry accepted: %s", accepted)
    return report


def assess_reference_hessian(directory: Path) -> dict[str, Any]:
    """Validate Hessian identity and finite unprojected frequency inventory.

    Frequency signs are retained for later constrained-curvature analysis.
    This function never claims an unconstrained minimum or force-field validity.

    Args:
        directory: Directory containing atom_map.json, hessian.fchk and
            frequency.log from the same constrained reference.

    Returns:
        Raw Hessian and frequency prerequisite checks with all parsed signs.

    Raises:
        ValueError: If atom metadata or the formatted checkpoint is malformed.
        FileNotFoundError: If a required input file is absent.
    """
    metadata = _read_json(directory / "atom_map.json")
    if metadata.get("model_scope") != "bonded_fit_surrogate_only":
        raise ValueError("Expected the declared capped bonded surrogate")
    if metadata.get("native_esp_resp_transfer_eligible") is not False:
        raise ValueError("Reference map must exclude native ESP/RESP transfer")
    atomic_numbers = _ordered_atomic_numbers(metadata)
    data = parse_gaussian_fchk(
        (directory / "hessian.fchk").read_text(encoding="utf-8"),
        expected_atomic_numbers=atomic_numbers,
        expected_charge=-1,
        expected_multiplicity=1,
    )
    text = (directory / "frequency.log").read_text(encoding="utf-8")
    frequencies = _read_frequencies(text)
    expected_modes = 3 * len(atomic_numbers) - 6
    normal = _ends_with_normal_termination(text) and "Error termination" not in text
    accepted = bool(
        normal
        and len(frequencies) == expected_modes
        and bool(np.isfinite(frequencies).all())
    )
    report = {
        "accepted": accepted,
        "scope": "raw constrained-reference Hessian prerequisite only",
        "forcefield_validated": False,
        "native_charge_transfer_eligible": False,
        "conditional_curvature_validation_pending": True,
        "unconstrained_minimum_claimed": False,
        "hessian_shape": list(data.hessian_eh_per_bohr2.shape),
        "hessian_units": "Hartree/bohr^2, not mass weighted",
        "unprojected_frequencies_cm_inverse": frequencies,
        "expected_unprojected_modes": expected_modes,
    }
    LOGGER.info("TI1 constrained reference Hessian accepted: %s", accepted)
    return report


def _read_json(path: Path) -> dict[str, Any]:
    """Load one JSON object and reject non-object roots explicitly."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path.name}")
    return payload


def _read_one_sdf(path: Path) -> Chem.Mol:
    """Read one valid explicit-hydrogen molecule from an SDF file."""
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Expected one valid explicit-hydrogen reference molecule")
    return molecules[0]


def _ordered_atomic_numbers(metadata: dict[str, Any]) -> tuple[int, ...]:
    """Validate atom-origin ordering and return the declared atomic numbers."""
    atom_map = metadata.get("atom_index_by_origin")
    records = metadata.get("atoms")
    if not isinstance(atom_map, dict) or not isinstance(records, list):
        raise ValueError("Reference map lacks ordered atom identity records")
    if len(atom_map) != 37 or len(records) != 37:
        raise ValueError("Reference map must identify all 37 atoms uniquely")
    if set(atom_map) != {record.get("origin_id") for record in records}:
        raise ValueError("Atom map origins differ from atom records")
    numbers = []
    for index, record in enumerate(records):
        origin = record.get("origin_id")
        if record.get("index") != index or atom_map.get(origin) != index:
            raise ValueError("Atom-origin indexes are duplicated or inconsistent")
        number = Chem.GetPeriodicTable().GetAtomicNumber(record.get("element", ""))
        if number <= 0:
            raise ValueError("Reference atom map contains an unknown element")
        numbers.append(number)
    return tuple(numbers)


def _validate_model_identity(
    molecule: Chem.Mol, metadata: dict[str, Any]
) -> dict[str, int]:
    """Validate atom map, charge, named graph, and stereochemical identity."""
    if metadata.get("model_scope") != "bonded_fit_surrogate_only":
        raise ValueError("Expected the declared capped bonded surrogate")
    if metadata.get("native_esp_resp_transfer_eligible") is not False:
        raise ValueError("Reference map must exclude native ESP/RESP transfer")
    if molecule.GetNumAtoms() != 37 or Chem.GetFormalCharge(molecule) != -1:
        raise ValueError("Expected 37-atom charge -1 reference")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    atom_map = metadata.get("atom_index_by_origin")
    records = metadata.get("atoms")
    _ordered_atomic_numbers(metadata)
    if not _REQUIRED_ORIGINS <= atom_map.keys():
        raise ValueError("Reference map lacks required reactive atom origins")
    seen_indices: set[int] = set()
    for record in records:
        origin = record["origin_id"]
        index = atom_map[origin]
        if not isinstance(index, int) or not 0 <= index < molecule.GetNumAtoms():
            raise ValueError(f"Invalid atom index for {origin}")
        if index in seen_indices or record.get("index") != index:
            raise ValueError("Atom-origin indexes are duplicated or inconsistent")
        seen_indices.add(index)
        atom = molecule.GetAtomWithIdx(index)
        if atom.GetSymbol() != record.get("element"):
            raise ValueError(f"Element mismatch for {origin}")
        if atom.GetFormalCharge() != record.get("formal_charge"):
            raise ValueError(f"Formal charge mismatch for {origin}")
    if seen_indices != set(range(molecule.GetNumAtoms())):
        raise ValueError("Atom map does not cover the complete molecule")

    origins = {
        origin: molecule.GetAtomWithIdx(index) for origin, index in atom_map.items()
    }
    for origin, expected in (("THR:CA", "S"), ("THR:CB", "R")):
        atom = origins[origin]
        if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != expected:
            raise ValueError(f"{origin} stereochemistry must remain {expected}")
    attacked_cip = metadata.get("attacked_carbon_cip")
    attacked = origins["TIL:C*"]
    if attacked_cip not in {"R", "S"} or not attacked.HasProp("_CIPCode"):
        raise ValueError("Attacked carbon must have explicit R/S stereochemistry")
    if attacked.GetProp("_CIPCode") != attacked_cip:
        raise ValueError("Attacked-carbon CIP differs from atom-map declaration")
    _validate_reactive_connectivity(molecule, origins)
    return {str(origin): int(index) for origin, index in atom_map.items()}


def _validate_reactive_connectivity(
    molecule: Chem.Mol, atoms: dict[str, Chem.Atom]
) -> None:
    """Require the declared covalent TI1 graph and internal proton inventory."""
    for first, second in (
        ("THR:OG1", "TIL:C*"),
        ("TIL:C*", "TIL:O"),
        ("TIL:C*", "TIL:N"),
        ("TIL:C*", "TIL:acylCH3"),
    ):
        bond = molecule.GetBondBetweenAtoms(
            atoms[first].GetIdx(), atoms[second].GetIdx()
        )
        if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
            raise ValueError(f"Expected single reactive bond {first}-{second}")
    if atoms["TIL:O"].GetFormalCharge() != -1:
        raise ValueError("TIL:O must retain the oxyanion charge")
    if any(neighbor.GetSymbol() == "H" for neighbor in atoms["TIL:O"].GetNeighbors()):
        raise ValueError("TIL:O must not be covalently protonated")
    if any(neighbor.GetSymbol() == "H" for neighbor in atoms["THR:OG1"].GetNeighbors()):
        raise ValueError("Thr OG1 must remain deprotonated in the bonded surrogate")
    if not any(
        neighbor.GetSymbol() == "H" for neighbor in atoms["TIL:N"].GetNeighbors()
    ):
        raise ValueError("TIL:N must retain its explicit hydrogen")


def _validate_reference_contract(
    metadata: dict[str, Any], atom_map: dict[str, int]
) -> None:
    """Check source identity and exact declared cap/backbone torsion targets."""
    reference = metadata.get("coordinate_reference")
    if not isinstance(reference, dict):
        raise ValueError("Reference map lacks coordinate-source identity")
    if any(reference.get(key) != value for key, value in _EXPECTED_REFERENCE.items()):
        raise ValueError("Expected source reference chain A, Thr residue 233")
    if not np.isclose(
        float(reference.get("psi_degrees", np.nan)),
        _EXPECTED_PSI,
        atol=1e-6,
        rtol=0.0,
    ):
        raise ValueError("Source psi differs from the declared Thr A:233 reference")
    if not np.isclose(
        float(reference.get("chi1_degrees", np.nan)),
        _EXPECTED_CHI1,
        atol=1e-6,
        rtol=0.0,
    ):
        raise ValueError("Source chi1 differs from the declared Thr A:233 reference")
    cap_phi = metadata.get("cap_phi_override_degrees")
    if cap_phi not in _CAP_PHI_CHOICES:
        raise ValueError("Explicit cap phi must be -120 or -60 degrees")

    constraints = metadata.get("constraints")
    if not isinstance(constraints, list) or len(constraints) != 3:
        raise ValueError("Expected frozen phi, psi and chi1 constraints")
    by_name = {row.get("name"): row for row in constraints if isinstance(row, dict)}
    if set(by_name) != set(TORSION_LABELS):
        raise ValueError("Frozen constraints must define phi, psi and chi1 once")
    expected_targets = {
        "phi": float(cap_phi),
        "psi": _EXPECTED_PSI,
        "chi1": _EXPECTED_CHI1,
    }
    for name, labels in TORSION_LABELS.items():
        row = by_name[name]
        expected_atoms = [atom_map[label] + 1 for label in labels]
        if row.get("atoms_one_based") != expected_atoms:
            raise ValueError(f"Frozen {name} constraint atom order is incorrect")
        if not np.isclose(
            float(row.get("degrees", np.nan)),
            expected_targets[name],
            atol=1e-6,
            rtol=0.0,
        ):
            raise ValueError(
                f"Frozen {name} target differs from the declared reference"
            )


def _measure_constraints(
    conformer: Chem.Conformer, constraints: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Measure final frozen torsions and their wrapped angular errors."""
    torsions = []
    for constraint in constraints:
        indices = [int(index) - 1 for index in constraint["atoms_one_based"]]
        observed = float(rdMolTransforms.GetDihedralDeg(conformer, *indices))
        target = float(constraint["degrees"])
        error = abs((observed - target + 180.0) % 360.0 - 180.0)
        torsions.append(
            {**constraint, "observed_degrees": observed, "error_degrees": error}
        )
    return torsions


def _reactive_bond_distances(
    atom_map: dict[str, int], coordinates: np.ndarray
) -> dict[str, float]:
    """Measure the four predeclared tetrahedral-carbon bond distances."""
    carbon = atom_map["TIL:C*"]
    return {
        name: float(np.linalg.norm(coordinates[carbon] - coordinates[index]))
        for name, index in atom_map.items()
        if name in _BOND_WINDOWS
    }


def _oxyanion_hydrogen_contacts(
    atoms: list[dict[str, Any]],
    atom_map: dict[str, int],
    coordinates: np.ndarray,
) -> list[dict[str, Any]]:
    """Report nearby hydrogens without treating hydrogen bonds as failures."""
    oxygen = coordinates[atom_map["TIL:O"]]
    contacts = []
    for atom in atoms:
        if atom.get("element") != "H":
            continue
        distance = float(np.linalg.norm(oxygen - coordinates[int(atom["index"])]))
        if distance < 2.5:
            contacts.append({"atom": atom["origin_id"], "distance_A": distance})
    return contacts


def _ends_with_normal_termination(text: str) -> bool:
    """Return whether a nonempty Gaussian log ends with normal termination."""
    lines = text.rstrip().splitlines()
    return bool(
        lines
        and re.fullmatch(
            r"Normal termination of Gaussian 16(?:\s+at\s+.+)?",
            lines[-1].strip(),
        )
    )


def _read_frequencies(text: str) -> list[float]:
    """Parse all Gaussian frequency values, preserving positive and negative signs."""
    values: list[float] = []
    for line in re.findall(r"Frequencies --([^\n]+)", text):
        try:
            values.extend(float(value) for value in line.split())
        except ValueError as error:
            raise ValueError("Malformed frequency value in Gaussian log") from error
    return values
