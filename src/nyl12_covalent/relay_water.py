"""Search relay-water seeds against a fully hydrogenated Nyl12 topology."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from parmed.amber import AmberParm


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1.0e-10:
        raise ValueError("Cannot construct an axis from coincident atoms")
    return vector / norm


def _frame(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    reference = np.eye(3)[int(np.argmin(np.abs(axis)))]
    first = _unit(np.cross(axis, reference))
    return first, np.cross(axis, first)


def _contact_limits(topology: AmberParm, indices: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """Set geometric screening limits, not force-field parameters."""
    oxygen_limits, hydrogen_limits = [], []
    for index in indices:
        atom = topology.atoms[index]
        element = atom.atomic_number
        if element == 1:
            polar = any(p.atomic_number in {7, 8, 16} for p in atom.bond_partners)
            oxygen_limits.append(1.55 if polar else 2.0)
            hydrogen_limits.append(1.4)
        elif element in {7, 8}:
            oxygen_limits.append(2.5)
            hydrogen_limits.append(1.55)
        elif element == 6:
            oxygen_limits.append(2.75)
            hydrogen_limits.append(2.1)
        elif element == 16:
            oxygen_limits.append(2.95)
            hydrogen_limits.append(2.3)
        else:
            raise ValueError(f"Unsupported local element: {element}")
    return np.asarray(oxygen_limits), np.asarray(hydrogen_limits)


def search_relay_water(
    topology: AmberParm,
    coordinates: np.ndarray,
    oxyanion_name: str,
    retain: int = 5,
) -> dict[str, object]:
    """Find one-water seeds against all local protein and ligand hydrogens.

    Args:
        topology: TI1 topology with Thr929, LIG1313, WAT1314.
        coordinates: Coordinates in exact topology order, in Angstrom.
        oxyanion_name: Substrate oxygen atom name.
        retain: Number of spatially distinct water seeds to retain.

    Returns:
        Search settings, candidate coordinates, and explicit steric gates.
        A failed gate must block refinement.
    """
    if len(topology.residues) != 1314 or retain < 1:
        raise ValueError("Unexpected TI1 inventory or candidate count")
    thr = {a.name: a for a in topology.residues[928].atoms}
    ligand = {a.name: a for a in topology.residues[1312].atoms}
    water = {a.name: a for a in topology.residues[1313].atoms}
    if set(water) != {"O", "H1", "H2"}:
        raise ValueError("Expected one explicitly hydrogenated relay water")
    nitrogen, og1 = coordinates[thr["N"].idx], coordinates[thr["OG1"].idx]
    oxyanion = coordinates[ligand[oxyanion_name].idx]
    n_hydrogens = coordinates[[thr[n].idx for n in ("H1", "H2", "H3")]]
    midpoint = (nitrogen + og1) / 2
    selected = [
        a.idx for a in topology.atoms if a.residue.idx != 1313
        and np.linalg.norm(coordinates[a.idx] - midpoint) < 9.0
    ]
    environment = coordinates[selected]
    o_limits, h_limits = _contact_limits(topology, selected)
    axis = _unit(og1 - nitrogen)
    separation = float(np.linalg.norm(og1 - nitrogen))
    first, second = _frame(axis)
    oxygen_seeds = []
    for n_radius in np.linspace(2.6, 3.4, 9):
        for og_radius in np.linspace(2.6, 3.4, 9):
            along = (n_radius**2 - og_radius**2 + separation**2) / (2 * separation)
            squared_radius = n_radius**2 - along**2
            if squared_radius <= 0:
                continue
            center = nitrogen + along * axis
            radius = math.sqrt(squared_radius)
            for angle in np.linspace(0, 2 * math.pi, 72, endpoint=False):
                oxygen = center + radius * (math.cos(angle) * first + math.sin(angle) * second)
                overlap = np.maximum(o_limits - np.linalg.norm(environment - oxygen, axis=1), 0)
                donor = float(np.linalg.norm(n_hydrogens - oxygen, axis=1).min())
                score = (
                    100 * float(overlap @ overlap)
                    + (n_radius - 2.85) ** 2 + (og_radius - 2.85) ** 2
                    + 0.5 * (donor - 1.85) ** 2
                )
                oxygen_seeds.append((score, oxygen, n_radius, og_radius))
    oxygen_seeds.sort(key=lambda candidate: candidate[0])
    if not oxygen_seeds:
        raise ValueError("No oxygen sites satisfy the relay-distance search")
    candidates = []
    bond, angle = 0.9572, math.radians(104.52)
    for o_score, oxygen, n_radius, og_radius in oxygen_seeds[:240]:
        for label, acceptor in (("OG1", og1), ("oxyanion", oxyanion)):
            direction = _unit(acceptor - oxygen)
            u, v = _frame(direction)
            h1 = oxygen + bond * direction
            for roll in np.linspace(0, 2 * math.pi, 36, endpoint=False):
                h2 = oxygen + bond * (
                    math.cos(angle) * direction
                    + math.sin(angle) * (math.cos(roll) * u + math.sin(roll) * v)
                )
                h_dist = np.linalg.norm(
                    np.array([h1, h2])[:, None, :] - environment[None, :, :], axis=2
                )
                h_overlap = np.maximum(h_limits - h_dist, 0)
                o_overlap = np.maximum(o_limits - np.linalg.norm(environment - oxygen, axis=1), 0)
                maximum = max(float(h_overlap.max()), float(o_overlap.max()))
                candidates.append({
                    "score": float(o_score + 100 * np.square(h_overlap).sum()),
                    "accepted": maximum <= 0.10,
                    "maximum_contact_shortfall_A": maximum,
                    "water_xyz_A": np.array([oxygen, h1, h2]).tolist(),
                    "hydrogen_1_donor_target": label,
                    "n_to_water_o_A": float(n_radius),
                    "og1_to_water_o_A": float(og_radius),
                    "oxyanion_to_water_o_A": float(np.linalg.norm(oxygen - oxyanion)),
                })
    candidates.sort(key=lambda c: (not c["accepted"], c["score"]))
    distinct = []
    for candidate in candidates:
        point = np.asarray(candidate["water_xyz_A"][0])
        if all(np.linalg.norm(point - np.asarray(p["water_xyz_A"][0])) >= 0.4 for p in distinct):
            distinct.append(candidate)
        if len(distinct) == retain:
            break
    return {
        "method": "sphere-intersection oxygen grid and hydrogen torsion scan",
        "oxygen_distance_range_A": [2.6, 3.4],
        "oxygen_distance_step_A": 0.1,
        "oxygen_angles": 72,
        "hydrogen_rolls": 36,
        "local_environment_atom_count": len(selected),
        "oxygen_seed_count": len(oxygen_seeds),
        "oriented_candidate_count": len(candidates),
        "accepted_oriented_candidate_count": sum(bool(c["accepted"]) for c in candidates),
        "maximum_allowed_contact_shortfall_A": 0.10,
        "water_atom_indices_zero_based": [water[n].idx for n in ("O", "H1", "H2")],
        "candidates": distinct,
        "interpretation": "Geometric seeds; a viable water is not evidence of a proton-transfer pathway.",
    }
