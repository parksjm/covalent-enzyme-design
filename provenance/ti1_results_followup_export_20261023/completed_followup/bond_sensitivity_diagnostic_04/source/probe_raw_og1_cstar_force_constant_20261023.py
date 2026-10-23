"""Probe the installed modified-Seminario result for one raw TI1 Hessian block."""

from __future__ import annotations

import hashlib
import inspect
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
from pymsmt.mcpb import gene_final_frcmod_file as fitting
from pymsmt.mol import constants
from nyl12_covalent.gaussian_fchk import read_gaussian_fchk

LOGGER = logging.getLogger(__name__)
ROOT = Path("/data/jpc/covalent-enzyme-design")
BRANCH = "attack_R_cap_phi_minus60"
BRANCH_ROOT = ROOT / "output_nyl12_ti1_parameters_20261022/pa6_local_bonded" / BRANCH
OUT = ROOT / "output_nyl12_parameter_validation_20261023/bond_sensitivity_diagnostic_04"
ORIGINS = ("THR:OG1", "TIL:C*")


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    """Read a JSON object and reject other root types."""
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def json_number(value: float) -> dict[str, Any]:
    """Represent a floating-point component without hiding nonfinite values."""
    finite = bool(np.isfinite(value))
    return {
        "finite": finite,
        "value": float(value) if finite else None,
        "repr": repr(value),
    }


def raw_result(result: tuple[Any, ...]) -> dict[str, Any]:
    """Preserve tuple values and complex components without altering signs."""
    force_constant = complex(result[1])
    payload: dict[str, Any] = {
        "tuple_length": len(result),
        "tuple_types": [type(item).__name__ for item in result],
        "raw_tuple_repr": repr(result),
        "distance_angstrom": json_number(float(result[0])),
        "force_constant_type": type(result[1]).__name__,
        "force_constant_repr": repr(result[1]),
        "force_constant_real": json_number(float(force_constant.real)),
        "force_constant_imaginary": json_number(float(force_constant.imag)),
        "force_constant_real_sign": (
            "positive"
            if force_constant.real > 0.0
            else "negative"
            if force_constant.real < 0.0
            else "zero"
        ),
        "force_constant_finite": bool(
            np.isfinite(force_constant.real) and np.isfinite(force_constant.imag)
        ),
    }
    if len(result) >= 3:
        deviation = complex(result[2])
        payload["standard_deviation_repr"] = repr(result[2])
        payload["standard_deviation_real"] = json_number(float(deviation.real))
        payload["standard_deviation_imaginary"] = json_number(float(deviation.imag))
    return payload


def main() -> None:
    """Report raw routine outputs for the designated OG1–C* Hessian block."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    report_path = OUT / "report.json"
    if report_path.exists():
        raise FileExistsError(report_path)
    paths = {
        "optimization_acceptance": BRANCH_ROOT
        / "optimization_01/scientific_acceptance.json",
        "frequency_acceptance": BRANCH_ROOT / "frequency_01/scientific_acceptance.json",
        "frequency_process": BRANCH_ROOT / "frequency_01/process_result.json",
        "input_identity": BRANCH_ROOT / "frequency_01/input_identity.json",
        "atom_map": BRANCH_ROOT / "frequency_01/atom_map.json",
        "fchk": BRANCH_ROOT / "frequency_01/hessian.fchk",
        "curvature": BRANCH_ROOT / "conditional_curvature_01/assessment.json",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    opt = read_json(paths["optimization_acceptance"])
    freq = read_json(paths["frequency_acceptance"])
    curvature = read_json(paths["curvature"])
    atom_map = read_json(paths["atom_map"])
    if opt.get("accepted") is not True or freq.get("accepted") is not True:
        raise ValueError("Upstream optimization/frequency acceptance is missing")
    if freq.get("checkpoint_geometry_match", {}).get("accepted") is not True:
        raise ValueError("Frequency checkpoint geometry was not accepted")
    if curvature.get("accepted") is not True:
        raise ValueError("Conditional-curvature assessment was not accepted")
    if atom_map.get("attacked_carbon_cip") != "R":
        raise ValueError("Unexpected attacked-carbon configuration")
    input_identity = read_json(paths["input_identity"])
    for filename, digest in input_identity.items():
        identity_path = paths["fchk"].parent / filename
        if not identity_path.is_file() or sha256(identity_path) != digest:
            raise ValueError(f"Frequency identity mismatch: {filename}")
    expected_fchk_hash = curvature["input_source_sha256"][str(paths["fchk"].resolve())]
    if sha256(paths["fchk"]) != expected_fchk_hash:
        raise ValueError("FCHK does not match accepted curvature source hash")
    atom_indices = atom_map["atom_index_by_origin"]
    i, j = (atom_indices[origin] for origin in ORIGINS)
    atom_elements = {"H": 1, "C": 6, "N": 7, "O": 8}
    expected_z = tuple(atom_elements[item["element"]] for item in atom_map["atoms"])
    fchk = read_gaussian_fchk(
        paths["fchk"],
        expected_atomic_numbers=expected_z,
        expected_charge=-1,
        expected_multiplicity=1,
    )
    coords = fchk.coordinates_bohr.reshape(-1)
    hessian = fchk.hessian_eh_per_bohr2
    forward = fitting.get_bond_fc_with_sem([], coords, hessian, i + 1, j + 1, 1.0, 0)
    reverse = fitting.get_bond_fc_with_sem([], coords, hessian, j + 1, i + 1, 1.0, 0)
    averaged = fitting.get_bond_fc_with_sem([], coords, hessian, i + 1, j + 1, 1.0, 1)
    block_ij = -hessian[3 * i : 3 * i + 3, 3 * j : 3 * j + 3]
    block_ji = -hessian[3 * j : 3 * j + 3, 3 * i : 3 * i + 3]
    sources = [
        Path(__file__),
        Path(inspect.getfile(fitting)),
        Path(inspect.getfile(constants)),
        ROOT / "src/nyl12_covalent/gaussian_fchk.py",
        ROOT
        / "output_nyl12_parameter_validation_20261023/bond_sensitivity_diagnostic_03/source/diagnose_ti1_bond_sensitivity_20261023.py",
    ]
    report = {
        "scope": "Raw API behavior probe only; no parameter fitting, replacement, or scientific acceptance",
        "branch": BRANCH,
        "bond": "OG1-C*",
        "atom_origins": list(ORIGINS),
        "one_based_atom_indices": [i + 1, j + 1],
        "frequency_scaling": 1.0,
        "routine_average_blocks_argument": 1,
        "routine_tuple_convention": "(distance_A, force_constant, standard_deviation) when average_blocks=1",
        "energy_convention": "K*(r-r0)^2",
        "constants": {
            "bohr_to_angstrom": float(constants.B_TO_A),
            "hartree_per_bohr2_to_kcal_per_mol_angstrom2": float(
                constants.HB2_TO_KCAL_MOL_A2
            ),
        },
        "raw_results": {
            "directed_H_ij": raw_result(forward),
            "directed_H_ji": raw_result(reverse),
            "average_blocks_1": raw_result(averaged),
        },
        "cross_block_eigenvalues": {
            "minus_H_ij": [repr(value) for value in np.linalg.eigvals(block_ij)],
            "minus_H_ji": [repr(value) for value in np.linalg.eigvals(block_ji)],
        },
        "input_sha256": {str(path.resolve()): sha256(path) for path in paths.values()},
        "source_sha256": {str(path.resolve()): sha256(path) for path in sources},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    LOGGER.info("Wrote raw modified-Seminario API probe; no force field was fitted")


if __name__ == "__main__":
    main()
