"""Rigid geometry and explicit contact screens for TI1 relay repair."""

from __future__ import annotations

import math
from collections import deque
from typing import Any

import numpy as np

from .relay_water import _contact_limits

CONTACT_LIMITS_A = {
    "H-H": 1.40,
    "H-C": 2.10,
    "H-polar": 1.55,
    "H-nonpolar": 2.00,
    "H-S": 2.30,
    "C-C": 2.80,
    "C-polar": 2.75,
    "C-S": 3.00,
    "polar-polar": 2.50,
    "polar-S": 2.95,
    "S-S": 3.20,
}


def unit(vector: np.ndarray) -> np.ndarray:
    """Return a normalized vector, rejecting a degenerate axis."""
    norm = float(np.linalg.norm(vector))
    if norm < 1.0e-10:
        raise ValueError("Cannot construct an axis from coincident atoms")
    return np.asarray(vector, dtype=float) / norm


def rotate_points(
    coordinates: np.ndarray,
    indices: list[int] | np.ndarray,
    origin: np.ndarray,
    axis: np.ndarray,
    degrees: float,
) -> np.ndarray:
    """Return a copied coordinate array with selected atoms rigidly rotated."""
    result = np.asarray(coordinates, dtype=float).copy()
    selected = np.asarray(indices, dtype=int)
    direction = unit(axis)
    angle = math.radians(degrees)
    vector = result[selected] - np.asarray(origin, dtype=float)
    result[selected] = (
        vector * math.cos(angle)
        + np.cross(direction, vector) * math.sin(angle)
        + np.outer(vector @ direction, direction) * (1.0 - math.cos(angle))
        + np.asarray(origin, dtype=float)
    )
    return result


def component_after_cut(atom: Any, cut_partner: Any) -> list[int]:
    """Return a graph component after cutting one undirected bond."""
    found: set[int] = set()
    queue: deque[Any] = deque([atom])
    while queue:
        current = queue.popleft()
        if current.idx in found:
            continue
        found.add(current.idx)
        for partner in current.bond_partners:
            if (current is atom and partner is cut_partner) or (
                current is cut_partner and partner is atom
            ):
                continue
            if partner.idx not in found:
                queue.append(partner)
    return sorted(found)


def graph_exclusions(atom: Any, maximum_graph_distance: int = 2) -> set[int]:
    """Return atom indices separated by one or two topology bonds."""
    distances = {atom.idx: 0}
    queue: deque[Any] = deque([atom])
    while queue:
        current = queue.popleft()
        distance = distances[current.idx]
        if distance == maximum_graph_distance:
            continue
        for partner in current.bond_partners:
            if partner.idx not in distances:
                distances[partner.idx] = distance + 1
                queue.append(partner)
    return set(distances)


def axis_local_environment_indices(
    coordinates: np.ndarray,
    environment_indices: list[int],
    origin: np.ndarray,
    axis: np.ndarray,
    radial_motion_bound_A: float,
    *,
    maximum_contact_A: float = 3.20,
) -> list[int]:
    """Keep every atom that could contact an atom within the motion bound.

    Any moved point remains within ``radial_motion_bound_A`` of the rotation
    axis. An environment atom farther than that bound plus the largest contact
    threshold cannot enter contact under any rotation around the axis.
    """
    indices = np.asarray(environment_indices, dtype=int)
    direction = unit(axis)
    relative = coordinates[indices] - np.asarray(origin, dtype=float)
    axial = relative @ direction
    perpendicular = relative - np.outer(axial, direction)
    distances = np.linalg.norm(perpendicular, axis=1)
    bound = radial_motion_bound_A + maximum_contact_A
    return indices[distances <= bound].tolist()


def _contact_class(element: int) -> str:
    if element == 6:
        return "C"
    if element in {7, 8}:
        return "polar"
    if element == 16:
        return "S"
    if element == 1:
        return "H"
    raise ValueError(f"Unsupported local element: {element}")


def _hydrogen_is_polar(atom: Any) -> bool:
    return any(partner.atomic_number in {7, 8, 16} for partner in atom.bond_partners)


def pair_contact_limit(first: Any, second: Any) -> float:
    """Return the requested atom-pair geometric screening threshold."""
    a, b = (
        _contact_class(int(first.atomic_number)),
        _contact_class(int(second.atomic_number)),
    )
    if a == "H" and b == "H":
        key = "H-H"
    elif a == "H" or b == "H":
        hydrogen, other = (first, second) if a == "H" else (second, first)
        other_class = b if a == "H" else a
        if other_class == "C":
            key = "H-C"
        elif other_class == "S":
            key = "H-S"
        else:
            key = "H-polar" if _hydrogen_is_polar(hydrogen) else "H-nonpolar"
    elif a == "S" and b == "S":
        key = "S-S"
    elif a == "S" or b == "S":
        other = b if a == "S" else a
        key = f"{other}-S" if other == "C" else "polar-S"
    elif a == "C" and b == "C":
        key = "C-C"
    elif a == "C" or b == "C":
        key = "C-polar"
    else:
        key = "polar-polar"
    return CONTACT_LIMITS_A[key]


def moved_atom_contact_screen(
    topology: Any,
    coordinates: np.ndarray,
    moved_indices: list[int],
    environment_indices: list[int],
    *,
    tolerance: float = 0.10,
) -> dict[str, Any]:
    """Measure moved-atom clashes, omitting only 1-2 and 1-3 graph pairs."""
    maximum = 0.0
    violations: list[dict[str, Any]] = []
    for moved_idx in moved_indices:
        moved = topology.atoms[moved_idx]
        excluded = graph_exclusions(moved)
        for other_idx in environment_indices:
            if other_idx == moved_idx or other_idx in excluded:
                continue
            other = topology.atoms[other_idx]
            distance = float(
                np.linalg.norm(coordinates[moved_idx] - coordinates[other_idx])
            )
            limit = pair_contact_limit(moved, other)
            shortfall = max(0.0, limit - distance)
            maximum = max(maximum, shortfall)
            if shortfall > tolerance:
                violations.append(
                    {
                        "moved_atom_index_zero_based": moved_idx,
                        "environment_atom_index_zero_based": other_idx,
                        "distance_A": distance,
                        "threshold_A": limit,
                        "shortfall_A": shortfall,
                    }
                )
    return {
        "accepted": maximum <= tolerance,
        "maximum_shortfall_A": maximum,
        "violation_count": len(violations),
        "violations": violations,
    }


def oxygen_grid(
    nitrogen: np.ndarray, og1: np.ndarray
) -> list[tuple[np.ndarray, float, float, int, int, int]]:
    """Generate all realizable points on the exact 9x9x72 distance grid."""
    axis_vector = np.asarray(og1, dtype=float) - np.asarray(nitrogen, dtype=float)
    separation = float(np.linalg.norm(axis_vector))
    axis = unit(axis_vector)
    reference = np.eye(3)[int(np.argmin(np.abs(axis)))]
    first = unit(np.cross(axis, reference))
    second = np.cross(axis, first)
    sites = []
    for i, n_distance in enumerate(np.linspace(2.6, 3.4, 9)):
        for j, og_distance in enumerate(np.linspace(2.6, 3.4, 9)):
            along = (n_distance**2 - og_distance**2 + separation**2) / (2 * separation)
            radius_squared = n_distance**2 - along**2
            if radius_squared <= 0:
                continue
            center = nitrogen + along * axis
            radius = math.sqrt(radius_squared)
            for k, angle in enumerate(
                np.linspace(0.0, 2 * math.pi, 72, endpoint=False)
            ):
                point = center + radius * (
                    math.cos(angle) * first + math.sin(angle) * second
                )
                sites.append((point, float(n_distance), float(og_distance), i, j, k))
    return sites


def contact_limits_for_indices(
    topology: Any, indices: list[int]
) -> tuple[np.ndarray, np.ndarray]:
    """Expose existing relay-water oxygen and hydrogen screen limits."""
    return _contact_limits(topology, indices)


def water_geometry(
    oxygen: np.ndarray, og1: np.ndarray, roll_degrees: float
) -> np.ndarray:
    """Build rigid water with H1 directed toward OG1 and a rolled H2."""
    direction = unit(np.asarray(og1, dtype=float) - oxygen)
    reference = np.eye(3)[int(np.argmin(np.abs(direction)))]
    first = unit(np.cross(direction, reference))
    second = np.cross(direction, first)
    angle = math.radians(104.52)
    roll = math.radians(roll_degrees)
    h1 = oxygen + 0.9572 * direction
    h2 = oxygen + 0.9572 * (
        math.cos(angle) * direction
        + math.sin(angle) * (math.cos(roll) * first + math.sin(roll) * second)
    )
    return np.asarray([oxygen, h1, h2])


def donor_geometry(
    nitrogen: np.ndarray, hydrogen: np.ndarray, acceptor: np.ndarray
) -> tuple[float, float]:
    """Return H···acceptor distance and N-H···acceptor angle in degrees."""
    h_to_n = unit(np.asarray(nitrogen, dtype=float) - hydrogen)
    h_to_acceptor = unit(np.asarray(acceptor, dtype=float) - hydrogen)
    cosine = float(np.clip(h_to_n @ h_to_acceptor, -1.0, 1.0))
    return (
        float(np.linalg.norm(np.asarray(acceptor, dtype=float) - hydrogen)),
        math.degrees(math.acos(cosine)),
    )
