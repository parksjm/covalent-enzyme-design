"""Assess standalone Gaussian HF/MK ESP data without fitting charges."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem

from nyl12_covalent.gaussian_checks import final_orientation, validate_geometry

BOHR_ANGSTROM = 0.529177210903


def _read_espgen(path: Path, atom_count: int) -> tuple[NDArray, NDArray]:
    """Read one complete, finite, standard espgen coordinate/potential block."""
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError("Empty espgen data")
    header = lines[0].split()
    if len(header) != 3:
        raise ValueError("Expected espgen atom/point/zero header")
    atoms, points, extra = map(int, header)
    if atoms != atom_count or points < 1 or extra != 0:
        raise ValueError("Invalid espgen atom count, point count or extra block")
    if len(lines) != 1 + atoms + points:
        raise ValueError("Truncated or extra espgen records")
    rows = [line.replace("D", "E").replace("d", "e").split() for line in lines[1:]]
    if any(len(row) != 3 for row in rows[:atoms]) or any(
        len(row) != 4 for row in rows[atoms:]
    ):
        raise ValueError("Invalid espgen coordinate/potential column count")
    centers = np.asarray(rows[:atoms], dtype=float)
    potentials = np.asarray(rows[atoms:], dtype=float)
    if not np.isfinite(centers).all() or not np.isfinite(potentials).all():
        raise ValueError("Nonfinite espgen coordinate or potential")
    return centers, potentials


def assess_hf_esp(
    log_path: Path,
    molecule_path: Path,
    esp_path: Path,
    reference_coordinates_angstrom: NDArray[np.float64],
    *,
    multiplicity: int,
    coordinate_tolerance_angstrom: float,
) -> dict[str, Any]:
    """Check standalone restricted-HF/6-31G(d) MK data at a fixed geometry.

    The caller separately verifies process completion and input/source hashes.
    This checks the log, graph, electronic state, unchanged atom coordinates and
    all extracted espgen records; it never certifies a RESP fit or force field.

    Args:
        log_path: One standalone Gaussian HF/MK calculation log.
        molecule_path: One explicit-hydrogen SDF with the intended atom order.
        esp_path: Unmodified Amber espgen output, in atomic units.
        reference_coordinates_angstrom: Accepted source geometry in atom order.
        multiplicity: Explicit expected spin multiplicity; this RHF stage is singlet.
        coordinate_tolerance_angstrom: Positive tolerance for printed coordinates.

    Returns:
        Separate scientific checks, geometry/grid summary and acceptance flag.

    Raises:
        ValueError: If inputs, geometry or ESP data are missing or malformed.
    """
    if multiplicity != 1:
        raise ValueError("This restricted-HF ESP stage requires singlet multiplicity")
    if (
        not np.isfinite(coordinate_tolerance_angstrom)
        or coordinate_tolerance_angstrom <= 0
    ):
        raise ValueError("Coordinate tolerance must be finite and positive")
    molecules = list(Chem.SDMolSupplier(str(molecule_path), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Expected one valid explicit-hydrogen SDF molecule")
    molecule = molecules[0]
    if any(
        atom.GetNumImplicitHs() or atom.GetNumExplicitHs()
        for atom in molecule.GetAtoms()
    ):
        raise ValueError("Every hydrogen must be a separate graph atom")
    reference = np.asarray(reference_coordinates_angstrom, dtype=float)
    if (
        reference.shape != (molecule.GetNumAtoms(), 3)
        or not np.isfinite(reference).all()
    ):
        raise ValueError("Invalid reference coordinates")
    text = log_path.read_text()
    elements, xyz = final_orientation(text)
    geometry = validate_geometry(molecule, elements, xyz)
    charges = [
        (int(charge), int(spin))
        for charge, spin in re.findall(
            r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", text
        )
    ]
    energies = re.findall(r"SCF Done:\s+E\(RHF\)\s*=\s*([-+\d.DEde]+)", text)
    energy = (
        float(energies[-1].replace("D", "E").replace("d", "e")) if energies else None
    )
    centers, potentials = _read_espgen(esp_path, molecule.GetNumAtoms())
    log_delta = float(np.max(np.abs(xyz - reference)))
    grid_delta = float(np.max(np.abs(centers * BOHR_ANGSTROM - xyz)))
    compact = re.sub(r"\s+", "", text).lower()
    checks = {
        "standalone_normal_termination": (
            text.count("Normal termination of Gaussian 16") == 1
            and "Normal termination of Gaussian 16" in text.rstrip().splitlines()[-1]
            and "Error termination" not in text
            and "Optimization completed." not in text
        ),
        "method_and_basis": bool(re.search(r"#p\s+HF/6-31G\(d\)", text, re.I))
        and "Standard basis: 6-31G(d)" in text,
        "mk_output_requested": all(
            item in compact for item in ("pop=mk", "6/33=2", "6/42=6")
        ),
        "electronic_state": bool(charges)
        and all(
            item == (Chem.GetFormalCharge(molecule), multiplicity) for item in charges
        ),
        "finite_restricted_hf_energy": energy is not None and bool(np.isfinite(energy)),
        "esp_log_present": "ESP Fit Center" in text
        and "Electrostatic Properties" in text,
        "geometry": bool(geometry["accepted"]),
        "source_geometry_unchanged": log_delta <= coordinate_tolerance_angstrom,
        "esp_atom_centers_match_log": grid_delta <= coordinate_tolerance_angstrom,
    }
    return {
        "accepted": all(checks.values()),
        "scope": "Standalone capped-surrogate HF/MK ESP prerequisites only",
        "checks": checks,
        "forcefield_validated": False,
        "native_charge_transfer_eligible": False,
        "charges_fitted": False,
        "charge": Chem.GetFormalCharge(molecule),
        "multiplicity": multiplicity,
        "energy_hartree": energy,
        "atom_count": molecule.GetNumAtoms(),
        "esp_grid_points": len(potentials),
        "max_source_coordinate_delta_angstrom": log_delta,
        "max_esp_atom_center_delta_angstrom": grid_delta,
        "coordinate_tolerance_angstrom": coordinate_tolerance_angstrom,
        "geometry": geometry,
    }
