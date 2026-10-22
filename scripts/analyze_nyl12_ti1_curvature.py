"""Assess constrained TI1 reference stability from a completed real Hessian."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.constrained_curvature import reduced_lagrangian_curvature  # noqa: E402
from nyl12_covalent.dihedral_constraints import (  # noqa: E402
    dihedral_constraint_derivatives,
    rigid_cartesian_directions,
)
from nyl12_covalent.gaussian_checks import validate_geometry  # noqa: E402
from nyl12_covalent.gaussian_fchk import read_gaussian_fchk  # noqa: E402
from nyl12_covalent.ti1_reference_coordinates import TORSION_LABELS  # noqa: E402

ANGSTROM_PER_BOHR = 0.529177210903
LOGGER = logging.getLogger(__name__)


def analyze_reference_curvature(
    frequency_dir: Path,
    output: Path,
    *,
    step_bohr: float,
    stationarity_tolerance: float,
    curvature_tolerance: float,
    derivative_relative_tolerance: float,
) -> dict[str, Any]:
    """Assess an accepted capped-reference Hessian without fitting parameters.

    Args:
        frequency_dir: Completed frequency_01 directory from the gated workflow.
        output: New analysis directory; must not exist.
        step_bohr: Finite-difference constraint derivative step in bohr.
        stationarity_tolerance: Maximum gradient residual in Hartree/bohr.
        curvature_tolerance: Unweighted eigenvalue tolerance in Hartree/bohr^2.
        derivative_relative_tolerance: Maximum relative derivative change on halving step.

    Returns:
        Stability and numerical checks for the declared frozen torsions only.

    Raises:
        ValueError: If prerequisites, identity, geometry or numerical inputs fail.
        FileExistsError: If the output already exists.
    """
    if output.exists():
        raise FileExistsError(output)
    if (
        not np.isfinite(derivative_relative_tolerance)
        or derivative_relative_tolerance <= 0
    ):
        raise ValueError("Derivative tolerance must be finite and positive")
    paths = [
        frequency_dir / name
        for name in (
            "scientific_acceptance.json",
            "process_result.json",
            "input_identity.json",
            "atom_map.json",
            "molecule.sdf",
            "hessian.fchk",
        )
    ]
    raw = json.loads(paths[0].read_text())
    process = json.loads(paths[1].read_text())
    if (
        raw["accepted"] is not True
        or raw["checkpoint_geometry_match"]["accepted"] is not True
        or process["exit_status"] != 0
    ):
        raise ValueError(
            "Completed and accepted raw Hessian/geometry prerequisites are required"
        )
    identities = json.loads(paths[2].read_text())
    required_inputs = {"reference.chk", "atom_map.json", "molecule.sdf", "input.gjf"}
    if not isinstance(identities, dict) or not required_inputs <= identities.keys():
        raise ValueError("Frequency input identity is incomplete")
    for filename, expected in identities.items():
        with (frequency_dir / filename).open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if actual != expected:
            raise ValueError(f"Frequency input identity changed: {filename}")
    metadata = json.loads(paths[3].read_text())
    if (
        metadata["model_scope"] != "bonded_fit_surrogate_only"
        or metadata["geometry_method"] != "wB97XD/def2TZVP"
    ):
        raise ValueError("Wrong model scope or geometry method")
    if metadata["native_esp_resp_transfer_eligible"] is not False:
        raise ValueError("Native charge transfer must remain ineligible")
    molecules = list(Chem.SDMolSupplier(str(paths[4]), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None or molecules[0].GetNumAtoms() != 37:
        raise ValueError("Expected one 37-atom molecular graph")
    molecule = molecules[0]
    data = read_gaussian_fchk(
        paths[5],
        expected_atomic_numbers=tuple(a.GetAtomicNum() for a in molecule.GetAtoms()),
        expected_charge=-1,
        expected_multiplicity=1,
    )
    if data.gradient_eh_per_bohr is None:
        raise ValueError(
            "Cartesian gradient is required; no zero-gradient substitution"
        )
    xyz_angstrom = data.coordinates_bohr * ANGSTROM_PER_BOHR
    geometry = validate_geometry(molecule, data.atomic_numbers.tolist(), xyz_angstrom)
    if not geometry["accepted"]:
        raise ValueError("Hessian coordinates fail molecular geometry screening")
    constraints = metadata["constraints"]
    if len(constraints) != 3 or {row["name"] for row in constraints} != {
        "phi",
        "psi",
        "chi1",
    }:
        raise ValueError("Expected the three declared frozen torsions")
    conformer = Chem.Conformer(37)
    for index, point in enumerate(xyz_angstrom):
        conformer.SetAtomPosition(index, tuple(float(x) for x in point))
    torsions = [
        tuple(index - 1 for index in row["atoms_one_based"]) for row in constraints
    ]
    for row, atoms in zip(constraints, torsions, strict=True):
        expected_atoms = tuple(
            metadata["atom_index_by_origin"][label]
            for label in TORSION_LABELS[row["name"]]
        )
        if atoms != expected_atoms:
            raise ValueError("Frozen torsion atom mapping differs from reference")
        observed = rdMolTransforms.GetDihedralDeg(conformer, *atoms)
        if abs((observed - row["degrees"] + 180.0) % 360.0 - 180.0) > 0.10:
            raise ValueError("Hessian coordinates violate a frozen torsion")
    coarse_j, coarse_h = dihedral_constraint_derivatives(
        data.coordinates_bohr, torsions, step_bohr=step_bohr
    )
    fine_j, fine_h = dihedral_constraint_derivatives(
        data.coordinates_bohr, torsions, step_bohr=step_bohr / 2
    )
    relative_errors = {}
    for name, coarse, fine in (
        ("jacobian", coarse_j, fine_j),
        ("constraint_hessians", coarse_h, fine_h),
    ):
        scale = float(np.linalg.norm(fine))
        if scale == 0:
            raise ValueError("Zero norm for nonlinear constraint derivatives")
        relative_errors[name] = float(np.linalg.norm(coarse - fine) / scale)
    rigid = rigid_cartesian_directions(data.coordinates_bohr)
    assessments = [
        reduced_lagrangian_curvature(
            data.hessian_eh_per_bohr2,
            data.gradient_eh_per_bohr.ravel(),
            jacobian,
            second,
            excluded_directions=rigid,
            stationarity_tolerance=stationarity_tolerance,
            curvature_tolerance=curvature_tolerance,
        )
        for jacobian, second in ((coarse_j, coarse_h), (fine_j, fine_h))
    ]
    derivative_pass = max(relative_errors.values()) <= derivative_relative_tolerance
    report = {
        "accepted": bool(
            derivative_pass and all(row["accepted"] for row in assessments)
        ),
        "scope": "Constrained reference stability only; no parameter or free-minimum certification",
        "forcefield_validated": False,
        "native_charge_transfer_eligible": False,
        "unconstrained_minimum_claimed": False,
        "constraint_step_bohr": step_bohr,
        "derivative_relative_errors": relative_errors,
        "derivative_relative_tolerance": derivative_relative_tolerance,
        "derivative_step_check_passed": derivative_pass,
        "coarse": assessments[0],
        "fine": assessments[1],
    }
    source_paths = paths + [
        Path(__file__),
        ROOT / "src/nyl12_covalent/constrained_curvature.py",
        ROOT / "src/nyl12_covalent/dihedral_constraints.py",
        ROOT / "src/nyl12_covalent/gaussian_fchk.py",
    ]
    report["input_source_sha256"] = {
        str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in source_paths
    }
    output.mkdir(parents=True)
    (output / "assessment.json").write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info("Constrained reference stability accepted: %s", report["accepted"])
    return report


def main() -> None:
    """Parse explicit reference and numerical settings, then assess stability."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frequency-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--step-bohr", required=True, type=float)
    parser.add_argument("--stationarity-tolerance", required=True, type=float)
    parser.add_argument("--curvature-tolerance", required=True, type=float)
    parser.add_argument("--derivative-relative-tolerance", required=True, type=float)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    report = analyze_reference_curvature(
        args.frequency_dir,
        args.output,
        step_bohr=args.step_bohr,
        stationarity_tolerance=args.stationarity_tolerance,
        curvature_tolerance=args.curvature_tolerance,
        derivative_relative_tolerance=args.derivative_relative_tolerance,
    )
    if not report["accepted"]:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
