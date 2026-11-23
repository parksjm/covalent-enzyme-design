"""Geometry measurements for Nyl12 catalytic-site alignments."""

from __future__ import annotations

from typing import Any

import numpy as np

CONTACTS = (
    ("Thr233.OG1", "Asp272.OD1"),
    ("Thr233.OG1", "Asp272.OD2"),
    ("Thr233.OG1", "Asn186.ND2"),
    ("Phe113.N", "Asp272.OD1"),
    ("Phe113.N", "Asp272.OD2"),
    ("Asn186.ND2", "Asp272.OD1"),
    ("Asn186.ND2", "Asp272.OD2"),
    ("Asp272.OD1", "Asp274.OD1"),
    ("Asp272.OD1", "Asp274.OD2"),
    ("Asp272.OD2", "Asp274.OD1"),
    ("Asp272.OD2", "Asp274.OD2"),
)

AtomMap = dict[tuple[int, str], tuple[str, np.ndarray]]


def fit_coordinates(
    mobile: np.ndarray,
    target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return row-vector Kabsch rotation, centers, and RMSD."""
    if (
        mobile.shape != target.shape
        or mobile.ndim != 2
        or mobile.shape[1] != 3
        or mobile.shape[0] < 3
    ):
        raise ValueError("Alignment requires at least three paired 3D coordinates")
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    left, _, right = np.linalg.svd(
        (mobile - mobile_center).T @ (target - target_center)
    )
    rotation = left @ right
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right
    aligned = (mobile - mobile_center) @ rotation + target_center
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - target) ** 2, axis=1))))
    return rotation, mobile_center, target_center, rmsd


def _dihedral(
    first: np.ndarray,
    second: np.ndarray,
    third: np.ndarray,
    fourth: np.ndarray,
) -> float:
    """Compute a signed dihedral angle in degrees."""
    first_bond = first - second
    central_bond = third - second
    final_bond = fourth - third
    norm = np.linalg.norm(central_bond)
    if norm < 1e-12:
        raise ValueError("Degenerate dihedral central bond")
    axis = central_bond / norm
    first_plane = first_bond - np.dot(first_bond, axis) * axis
    final_plane = final_bond - np.dot(final_bond, axis) * axis
    if np.linalg.norm(first_plane) < 1e-12 or np.linalg.norm(final_plane) < 1e-12:
        raise ValueError("Degenerate dihedral plane")
    return float(
        np.degrees(
            np.arctan2(
                np.dot(np.cross(axis, first_plane), final_plane),
                np.dot(first_plane, final_plane),
            )
        )
    )


def wrapped_delta(first: float, second: float) -> float:
    """Return the signed wrapped difference ``first - second`` in degrees."""
    return float((first - second + 180.0) % 360.0 - 180.0)


def _point(source: AtomMap, label: str) -> np.ndarray:
    """Get one named catalytic atom coordinate."""
    residue, atom = label.split(".")
    residue_number = {
        "Phe113": 113,
        "Asn186": 186,
        "Thr233": 233,
        "Thr234": 234,
        "Asp272": 272,
        "Asp274": 274,
    }[residue]
    key = (residue_number, atom)
    if key not in source:
        raise ValueError(f"Missing contact atom {label}")
    return source[key][1]


def state_metrics(reference: AtomMap, model: AtomMap) -> dict[str, Any]:
    """Measure held-out Thr233 and active-site contacts for one aligned state."""
    thr_displacements = {
        atom: float(np.linalg.norm(model[(233, atom)][1] - reference[(233, atom)][1]))
        for auth, atom in model
        if auth == 233 and (233, atom) in reference
    }
    if not thr_displacements:
        raise ValueError("No matching heavy atoms found for held-out Thr233")
    contacts = {}
    for first, second in CONTACTS:
        reference_distance = float(
            np.linalg.norm(_point(reference, first) - _point(reference, second))
        )
        model_distance = float(
            np.linalg.norm(_point(model, first) - _point(model, second))
        )
        contacts[f"{first}--{second}"] = {
            "reference_angstrom": reference_distance,
            "aligned_es_angstrom": model_distance,
            "delta_es_minus_reference_angstrom": (model_distance - reference_distance),
        }
    chi_atoms = ("Thr233.N", "Thr233.CA", "Thr233.CB", "Thr233.OG1")
    psi_atoms = ("Thr233.N", "Thr233.CA", "Thr233.C", "Thr234.N")
    reference_chi = _dihedral(*(_point(reference, label) for label in chi_atoms))
    model_chi = _dihedral(*(_point(model, label) for label in chi_atoms))
    reference_psi = _dihedral(*(_point(reference, label) for label in psi_atoms))
    model_psi = _dihedral(*(_point(model, label) for label in psi_atoms))
    return {
        "thr233_held_out": {
            "heavy_atom_displacements_angstrom": thr_displacements,
            "rmsd_angstrom": float(
                np.sqrt(np.mean(np.square(list(thr_displacements.values()))))
            ),
            "maximum_displacement_angstrom": max(thr_displacements.values()),
            "chi1_reference_degrees": reference_chi,
            "chi1_aligned_es_degrees": model_chi,
            "chi1_delta_es_minus_reference_degrees": wrapped_delta(
                model_chi,
                reference_chi,
            ),
            "psi_reference_degrees": reference_psi,
            "psi_aligned_es_degrees": model_psi,
            "psi_delta_es_minus_reference_degrees": wrapped_delta(
                model_psi,
                reference_psi,
            ),
        },
        "contacts": contacts,
    }
