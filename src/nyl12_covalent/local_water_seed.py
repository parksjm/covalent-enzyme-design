"""Bounded geometric preparation of explicit water before energy refinement."""

from __future__ import annotations

from itertools import product

import numpy as np

from .relay_repair_geometry import unit


def seed_water(
    oxygen: np.ndarray,
    nitrogen: np.ndarray,
    og1: np.ndarray,
    environment: np.ndarray,
    elements: np.ndarray,
    *,
    preserve_oxygen: bool = False,
) -> tuple[np.ndarray, dict]:
    """Remove gross oxygen overlap and construct an oriented neutral H2O seed.

    Args:
        oxygen: Original AF3 oxygen coordinate.
        nitrogen: Catalytic terminal nitrogen coordinate.
        og1: Catalytic threonine OG1 coordinate.
        environment: All nonwater coordinates, including hydrogens.
        elements: Atomic numbers corresponding to environment coordinates.
        preserve_oxygen: Keep the original AF3 oxygen without geometric repair.

    Returns:
        O/H1/H2 coordinates and geometric provenance. This is not minimization.
    """
    if environment.shape != (len(elements), 3) or not np.isfinite(environment).all():
        raise ValueError("Invalid water environment")
    heavy = elements > 1
    thresholds = {6: 2.75, 7: 2.50, 8: 2.50, 16: 2.95}
    limits = np.array([thresholds[int(z)] for z in elements[heavy]])
    best = None
    count = 0
    offsets = (
        [(0.0, 0.0, 0.0)]
        if preserve_oxygen
        else product(np.arange(-1.5, 1.501, 0.25), repeat=3)
    )
    for offset in offsets:
        shift = np.asarray(offset)
        if np.linalg.norm(shift) > 1.5:
            continue
        trial = oxygen + shift
        bridge = [np.linalg.norm(trial - atom) for atom in (nitrogen, og1)]
        if not all(2.6 <= d <= 3.2 for d in bridge):
            continue
        distances = np.linalg.norm(environment[heavy] - trial, axis=1)
        shortfall = np.maximum(limits - distances, 0)
        score = (
            float(shortfall.max()),
            float(shortfall @ shortfall),
            float(shift @ shift),
        )
        count += 1
        if best is None or score < best[0]:
            best = (score, trial)
    if best is None:
        raise ValueError("No bounded water seed within both bridge ranges")
    score, oxygen_new = best
    minimum = float(np.linalg.norm(environment[heavy] - oxygen_new, axis=1).min())
    if minimum < (1e-6 if preserve_oxygen else 1.8):
        raise ValueError("Bounded seed retains a catastrophic heavy-atom overlap")
    axis = unit(og1 - oxygen_new)
    basis = np.eye(3)[np.argmin(np.abs(axis))]
    perpendicular = unit(np.cross(axis, basis))
    second = np.cross(axis, perpendicular)
    length, angle = 0.9572, np.deg2rad(104.52)
    candidates = []
    h_limits = np.where(elements == 1, 1.4, np.where(elements == 6, 2.1, 1.55))
    for degrees in range(0, 360, 30):
        phi = np.deg2rad(degrees)
        direction = np.cos(angle) * axis + np.sin(angle) * (
            np.cos(phi) * perpendicular + np.sin(phi) * second
        )
        water = np.array(
            [oxygen_new, oxygen_new + length * axis, oxygen_new + length * direction]
        )
        separations = np.linalg.norm(
            water[1:, None, :] - environment[None, :, :], axis=2
        )
        penalty = np.maximum(h_limits - separations, 0)
        candidates.append((float(np.sum(penalty**2)), degrees, water))
    penalty, degrees, water = min(candidates, key=lambda item: item[:2])
    return water, {
        "method": (
            "original AF3 oxygen unchanged; 12 H2 orientations; no energy evaluation"
            if preserve_oxygen
            else "bounded 0.25 A oxygen grid; 12 H2 orientations; no energy evaluation"
        ),
        "original_oxygen_A": oxygen.tolist(),
        "oxygen_displacement_A": float(np.linalg.norm(oxygen_new - oxygen)),
        "eligible_oxygen_candidates": count,
        "maximum_heavy_contact_shortfall_A": score[0],
        "minimum_heavy_distance_A": minimum,
        "hydrogen_orientation_degrees": degrees,
        "hydrogen_overlap_penalty_A2": penalty,
        "OH_length_A": length,
        "HOH_angle_degrees": 104.52,
        "accepted_minimum": False,
    }
