"""Geometry and input contracts for bounded N-terminal relay releases."""

from __future__ import annotations

from itertools import combinations
import re
from typing import Any

import numpy as np

from nyl12_covalent.native_sander_cli import prepare_native_input

ROLES = {
    "C": 19916,
    "O": 19917,
    "OG": 13916,
    "N": 13905,
    "NL": 19929,
    "H1": 13906,
    "H2": 13907,
    "HG": 13917,
    "HN": 19937,
    "CA": 19915,
}


def geometry(xyz: np.ndarray) -> dict[str, Any]:
    """Measure bonds and all relay proton alternatives without atom relabeling.

    Args:
        xyz: Full-system coordinates in original topology order, in Angstrom.

    Returns:
        Distances, carbonyl-center angles and a conservative geometry label.
    """
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < max(ROLES.values()):
        raise ValueError("Full original-topology coordinates required")
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite coordinates")
    pos = {key: xyz[value - 1] for key, value in ROLES.items()}
    pairs = [("C", atom) for atom in ("O", "OG", "NL", "CA")]
    pairs += [(h, a) for h in ("H1", "H2", "HG", "HN") for a in ("N", "OG", "NL", "O")]
    distances = {f"{a}_{b}": float(np.linalg.norm(pos[a] - pos[b])) for a, b in pairs}
    angles = {}
    for a, b in combinations(("O", "OG", "NL", "CA"), 2):
        u, v = pos[a] - pos["C"], pos[b] - pos["C"]
        denominator = np.linalg.norm(u) * np.linalg.norm(v)
        if denominator == 0:
            raise ValueError("Coincident carbon-center coordinates")
        angles[f"{a}_C_{b}"] = float(
            np.degrees(np.arccos(np.clip(np.dot(u, v) / denominator, -1, 1)))
        )
    if distances["C_OG"] < 1.7 and distances["C_NL"] < 1.7:
        label = "TI1_like_geometric_candidate"
    elif distances["C_OG"] > 2.3 and distances["C_NL"] < 1.7:
        label = "reactant_like"
    elif distances["C_OG"] < 1.7 and distances["C_NL"] > 2.0:
        label = "acyl_like"
    else:
        label = "unassigned"
    return {
        "distances_A": distances,
        "carbon_angles_degrees": angles,
        "geometry_label": label,
        "proton_bond_cutoff_A": 1.2,
        "proton_bonds": {
            h: [a for a in ("N", "OG", "NL", "O") if distances[f"{h}_{a}"] < 1.2]
            for h in ("H1", "H2", "HG", "HN")
        },
    }


def minimization_input(source: str, maxcycles: int = 10000) -> str:
    """Build native LBFGS input preserving the exact QM and movable masks.

    Args:
        source: Validated unbiased native dimer MDIN with immutable active mask.
        maxcycles: Positive maximum native minimizer iterations.

    Returns:
        Zero-temperature minimization input with no coordinate/position potential.
    """
    if maxcycles < 1:
        raise ValueError("Positive maxcycles required")
    for key in ("ntr", "nmropt"):
        found = re.findall(rf"\b{key}\s*=\s*(\d+)", source)
        if found != ["0"]:
            raise ValueError(f"Require exactly one {key}=0")
    result = prepare_native_input(source)
    result = re.sub(r"\bmaxcyc=1\b", f"maxcyc={maxcycles}", result)
    result = result.replace("maxcycle=1,", f"maxcycle={maxcycles},")
    return result.replace("tol=1.d-4, tole=1.d-5", "tol=1.d-5, tole=1.d-7")


def aligned_rmsd(left: np.ndarray, right: np.ndarray) -> float:
    """Return proper-rotation RMSD between equal finite coordinate selections."""
    if left.shape != right.shape or left.ndim != 2 or left.shape[1] != 3:
        raise ValueError("Equal N-by-3 selections required")
    if not np.isfinite(left).all() or not np.isfinite(right).all() or not len(left):
        raise ValueError("Finite nonempty coordinates required")
    a, b = left - left.mean(axis=0), right - right.mean(axis=0)
    u, _, vh = np.linalg.svd(a.T @ b)
    correction = np.diag([1.0, 1.0, np.linalg.det(u @ vh)])
    return float(np.sqrt(np.mean(np.sum((a @ u @ correction @ vh - b) ** 2, axis=1))))
