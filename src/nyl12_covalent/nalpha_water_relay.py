"""Exact periodic geometry for observed N-alpha/water/leaving-N relays."""

from __future__ import annotations

from itertools import product
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from nyl12_covalent.catalytic_alignment import BACKBONE, FIT_RESIDUES, SITE_RESIDUES
from nyl12_covalent.catalytic_alignment_geometry import fit_coordinates, state_metrics


def cell_matrix(lengths: np.ndarray, angles: np.ndarray) -> np.ndarray:
    """Return row-vector triclinic lattice from Amber lengths and angles."""
    a, b, c = np.asarray(lengths, dtype=float)
    alpha, beta, gamma = np.radians(angles)
    matrix = np.array(
        [
            [a, 0.0, 0.0],
            [b * np.cos(gamma), b * np.sin(gamma), 0.0],
            [
                c * np.cos(beta),
                c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
                0.0,
            ],
        ]
    )
    matrix[2, 2] = np.sqrt(c * c - np.dot(matrix[2], matrix[2]))
    if not np.isfinite(matrix).all() or np.linalg.det(matrix) <= 0:
        raise ValueError("Invalid periodic cell")
    return matrix


def minimum_image(vectors: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Return exact nearest lattice images using a singular-value bound.

    Fractional rounding supplies an upper bound only. Every integer translation
    capable of improving it is enumerated, including skew triclinic cells.
    """
    values = np.asarray(vectors, dtype=float)
    if values.shape[-1] != 3 or not np.isfinite(values).all():
        raise ValueError("Expected finite three-dimensional vectors")
    fractional = values.reshape(-1, 3) @ np.linalg.inv(cell)
    fractional -= np.rint(fractional)
    rounded = fractional @ cell
    if not len(rounded):
        return rounded.reshape(values.shape)
    radius = float(np.linalg.norm(rounded, axis=1).max())
    singular = float(np.linalg.svd(cell, compute_uv=False).min())
    bound = int(np.floor(radius / singular + 0.5))
    lattice = np.array(list(product(range(-bound, bound + 1), repeat=3))) @ cell
    # A KD tree avoids materializing vectors x lattice x xyz.
    nearest = cKDTree(lattice).query(rounded)[1]
    return (rounded - lattice[nearest]).reshape(values.shape)


def close_waters(
    positions: np.ndarray, anchor: np.ndarray, cell: np.ndarray, cutoff: float
) -> np.ndarray:
    """Return indices within cutoff, with an exact conservative fractional filter."""
    inverse = np.linalg.inv(cell)
    limits = cutoff * np.linalg.norm(inverse, axis=0)
    if np.any(limits >= 0.5):
        raise ValueError("Cell too small for the local-water fractional prefilter")
    fractional = (positions - anchor) @ inverse
    fractional -= np.rint(fractional)
    possible = np.flatnonzero(np.all(np.abs(fractional) <= limits, axis=1))
    displacement = minimum_image(fractional[possible] @ cell, cell)
    return possible[np.linalg.norm(displacement, axis=1) <= cutoff]


def angle(first: np.ndarray, vertex: np.ndarray, last: np.ndarray) -> float:
    """Return a nondegenerate three-point angle in degrees."""
    left, right = first - vertex, last - vertex
    norm = np.linalg.norm(left) * np.linalg.norm(right)
    if norm < 1e-12:
        raise ValueError("Degenerate angle")
    return float(np.degrees(np.arccos(np.clip(np.dot(left, right) / norm, -1, 1))))


def hbond(
    xyz: np.ndarray, donor: int, hydrogen: int, acceptor: int, cell: np.ndarray
) -> dict[str, float | bool]:
    """Measure a topology-owned D-H...A bond under exact periodic images."""
    dh, ha, da = minimum_image(
        np.array(
            [
                xyz[donor] - xyz[hydrogen],
                xyz[acceptor] - xyz[hydrogen],
                xyz[acceptor] - xyz[donor],
            ]
        ),
        cell,
    )
    distance = float(np.linalg.norm(ha))
    donor_distance = float(np.linalg.norm(da))
    degrees = angle(dh, np.zeros(3), ha)
    return {
        "pass": donor_distance <= 3.5 and distance <= 2.5 and degrees >= 120,
        "da_A": donor_distance,
        "ha_A": distance,
        "angle_deg": degrees,
    }


def native_geometry(
    xyz: np.ndarray, topology: Any, reference: dict, cell: np.ndarray, nalpha: int
) -> dict[str, float | bool]:
    """Apply existing held-out native Thr displacement and torsion gates."""
    indices, keys = [], []
    for amber, canonical, expected in SITE_RESIDUES:
        residue = topology.residues[amber - 1]
        allowed = {"ASP", "ASH"} if expected == "ASP" else {expected}
        if residue.name not in allowed:
            raise ValueError(f"Residue mapping mismatch: {amber} {residue.name}")
        for atom in residue.atoms:
            if atom.atomic_number != 1:
                indices.append(atom.idx)
                keys.append((canonical, atom.name))
    local = minimum_image(xyz[indices] - xyz[nalpha], cell) + xyz[nalpha]
    mobile = {
        key: (reference[key][0] if key in reference else "H", value)
        for key, value in zip(keys, local, strict=True)
    }
    fit_keys = [(res, atom) for _, res, _ in FIT_RESIDUES for atom in BACKBONE]
    rotation, center, target, rmsd = fit_coordinates(
        np.array([mobile[key][1] for key in fit_keys]),
        np.array([reference[key][1] for key in fit_keys]),
    )
    aligned = {
        key: (element, (point - center) @ rotation + target)
        for key, (element, point) in mobile.items()
    }
    thr = state_metrics(reference, aligned)["thr233_held_out"]
    shifts = thr["heavy_atom_displacements_angstrom"]
    psi = thr["psi_delta_es_minus_reference_degrees"]
    chi = thr["chi1_delta_es_minus_reference_degrees"]
    return {
        "native_thr_rmsd_A": thr["rmsd_angstrom"],
        "scaffold_rmsd_A": rmsd,
        "native_og1_shift_A": shifts["OG1"],
        "native_nalpha_shift_A": shifts["N"],
        "psi_delta_deg": psi,
        "chi1_delta_deg": chi,
        "native_pass": shifts["OG1"] <= 0.75
        and shifts["N"] <= 1.25
        and abs(psi) <= 35
        and abs(chi) <= 20,
    }


def relay_class(first: bool, second: bool) -> str:
    """Keep observed two-leg, leaving-leg and proximity-only classes separate."""
    if first and second:
        return "strict_both_leg_hbond"
    return "leaving_leg_hbond" if second else "heavy_bridge_only"


def candidate_rank(row: dict[str, Any]) -> tuple:
    """Rank observed chemistry classes before tier and continuous geometry."""
    classes = {
        "strict_both_leg_hbond": 0,
        "leaving_leg_hbond": 1,
        "heavy_bridge_only": 2,
    }
    return (
        classes[row["relay_class"]],
        0 if row["tier"] == "A" else 1,
        max(0, row["water_h_to_nleave_A"] - 2.5)
        + max(0, 120 - row["water_h_to_nleave_angle_deg"]) / 60,
        row["water_h_to_nleave_A"],
        row["native_thr_rmsd_A"],
        row["case"],
        row["frame"],
        row["water_resid"],
    )
