"""Build atom-preserving acylenzyme coordinate hypotheses for QM/MM only."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from parmed import Structure


@dataclass(frozen=True)
class AcylAtomMap:
    """One-based identities in the unchanged E-S topology."""

    thr_n: int = 13905
    thr_h1: int = 13906
    thr_h2: int = 13907
    thr_og: int = 13916
    proton: int = 13917
    carbon: int = 19916
    oxygen: int = 19917
    leaving_n: int = 19929


def _unit(vector: np.ndarray) -> np.ndarray:
    length = float(np.linalg.norm(vector))
    if length < 1e-8:
        raise ValueError("Degenerate coordinate direction")
    return vector / length


def _component(graph: list[set[int]], start: int) -> set[int]:
    visited, pending = set(), [start]
    while pending:
        atom = pending.pop()
        if atom not in visited:
            visited.add(atom)
            pending.extend(graph[atom] - visited)
    return visited


def _orient_leaving(
    coordinates: np.ndarray,
    moving: list[int],
    center: int,
    graph: list[set[int]],
    radii: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    fixed = np.array(sorted(set(range(len(coordinates))) - set(moving)))
    nearby = np.zeros(len(fixed), dtype=bool)
    for atom in moving:
        nearby |= np.linalg.norm(coordinates[fixed] - coordinates[atom], axis=1) < 6
    fixed = fixed[nearby]
    if not len(fixed):
        return coordinates, {"axis": None, "degrees": 0.0}
    valid = np.ones((len(moving), len(fixed)), dtype=bool)
    for row, atom in enumerate(moving):
        excluded = set(graph[atom])
        for neighbor in graph[atom]:
            excluded.update(graph[neighbor])
        valid[row] = [atom_index not in excluded for atom_index in fixed]
    cutoffs = 0.65 * (radii[moving, None] + radii[fixed][None, :])
    original = coordinates[moving] - coordinates[center]

    def score(points: np.ndarray) -> tuple[float, float]:
        distances = np.linalg.norm(points[:, None] - coordinates[fixed][None], axis=2)
        ratios = (distances / cutoffs)[valid]
        return (
            float(np.minimum(ratios, 1).min(initial=1)),
            -float(np.maximum(1 - ratios, 0).sum()),
        )

    best, best_points = score(coordinates[moving]), coordinates[moving].copy()
    selected: dict[str, Any] = {"axis": None, "degrees": 0.0}
    for axis_index, axis in enumerate(np.eye(3)):
        cross = np.array(
            [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
        )
        for degrees in (-5, 5, -10, 10, -15, 15):
            angle = np.radians(degrees)
            rotation = (
                np.eye(3)
                + np.sin(angle) * cross
                + (1 - np.cos(angle)) * (cross @ cross)
            )
            points = original @ rotation.T + coordinates[center]
            candidate = score(points)
            if candidate > best:
                best, best_points = candidate, points
                selected = {"axis": "xyz"[axis_index], "degrees": float(degrees)}
    result = coordinates.copy()
    result[moving] = best_points
    return result, selected


def _validated_graph(
    topology: Structure, mapping: AcylAtomMap
) -> tuple[list[set[int]], dict[str, int]]:
    names = ("N", "H1", "H2", "OG1", "HG1", "C", "O", "N")
    elements = (7, 1, 1, 8, 1, 6, 8, 7)
    indices = {name: value - 1 for name, value in vars(mapping).items()}
    if len(set(indices.values())) != len(indices):
        raise ValueError("Atom mapping must contain unique indices")
    for index, name, element in zip(indices.values(), names, elements, strict=True):
        if not 0 <= index < len(topology.atoms):
            raise ValueError("Atom mapping index outside topology")
        atom = topology.atoms[index]
        if atom.name != name or atom.atomic_number != element:
            raise ValueError(f"Atom mapping mismatch at one-based index {index + 1}")
    thr = [indices[name] for name in ("thr_n", "thr_h1", "thr_h2", "thr_og", "proton")]
    if len({topology.atoms[i].residue.idx for i in thr}) != 1:
        raise ValueError("Thr atom mapping spans residues")
    if topology.atoms[indices["thr_n"]].residue.name != "THR":
        raise ValueError("Mapped catalytic residue must be THR")
    graph = [{partner.idx for partner in atom.bond_partners} for atom in topology.atoms]
    required = (
        ("thr_n", "thr_h1"),
        ("thr_n", "thr_h2"),
        ("thr_og", "proton"),
        ("carbon", "oxygen"),
        ("carbon", "leaving_n"),
    )
    if any(indices[b] not in graph[indices[a]] for a, b in required):
        raise ValueError("Source E-S topology lacks required bond")
    if graph[indices["proton"]] != {indices["thr_og"]}:
        raise ValueError("HG1 must have exactly its source OG1 topology bond")
    return graph, indices


def build_acylenzyme_seed(
    topology: Structure,
    coordinates: np.ndarray,
    mapping: AcylAtomMap = AcylAtomMap(),
    *,
    ester_distance: float = 1.4,
    leaving_distance: float = 2.85,
    orient_leaving: bool = False,
    source_proton_state: Literal["thr_n", "leaving_n"] = "thr_n",
    leaving_bend_degrees: float = 0.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Construct a planar ester and a neutral, detached leaving amine.

    Args:
        topology: Unmodified E-S graph, atom order and identities.
        coordinates: Seed coordinates in angstrom, shape (atoms, 3).
        mapping: Explicit one-based reactive atom identities.
        ester_distance: Explicit initial OG1-C distance in angstrom.
        leaving_distance: Explicit initial broken C-N distance in angstrom.
        orient_leaving: Search deterministic rigid rotations up to 15 degrees
            about leaving N to reduce close pocket contacts, including HG1.
        source_proton_state: Explicit HG1 parent in source coordinates: Thr N
            for TI1 or leaving N for a direct transferred-proton candidate.
        leaving_bend_degrees: Explicit C-N direction rotation within the
            carbonyl plane, bounded to 30 degrees; fragment remains rigid.

    Returns:
        New coordinates and a JSON-compatible audit containing moved indices,
        fragment membership, reactive distances and close contacts. The source
        topology remains an E-S topology and is unsuitable for classical MD of
        the product. No energy minimization or scientific acceptance is implied.

    Raises:
        ValueError: If mapping, graph, source proton state or geometry is invalid.
    """
    source = np.asarray(coordinates, dtype=float)
    if source.shape != (len(topology.atoms), 3) or not np.isfinite(source).all():
        raise ValueError("Expected finite coordinates with shape (atoms, 3)")
    if not (1.3 <= ester_distance <= 1.8 and 2.3 <= leaving_distance <= 3.1):
        raise ValueError("Seed distances outside declared ester/cleaved-bond range")
    if not np.isfinite(leaving_bend_degrees) or abs(leaving_bend_degrees) > 30:
        raise ValueError("Leaving bend must be finite and within 30 degrees")
    graph, ids = _validated_graph(topology, mapping)
    n, og, h = (ids[key] for key in ("thr_n", "thr_og", "proton"))
    c, o, leaving = (ids[key] for key in ("carbon", "oxygen", "leaving_n"))
    if source_proton_state not in ("thr_n", "leaving_n"):
        raise ValueError("Unknown source proton state")
    for hydrogen in (ids["thr_h1"], ids["thr_h2"]):
        if not 0.8 <= np.linalg.norm(source[n] - source[hydrogen]) <= 1.3:
            raise ValueError("Source requires H1/H2 bound to Thr N")
    parent = n if source_proton_state == "thr_n" else leaving
    if not 0.8 <= np.linalg.norm(source[parent] - source[h]) <= 1.3:
        raise ValueError("HG1 geometry disagrees with explicit source proton state")
    if (
        source_proton_state == "leaving_n"
        and np.linalg.norm(source[n] - source[h]) <= 1.3
    ):
        raise ValueError("Transferred HG1 remains within bonding distance of Thr N")
    acyl_neighbors = graph[c] - {o, leaving}
    if len(acyl_neighbors) != 1 or graph[o] != {c}:
        raise ValueError("Expected carbonyl C bonded only to O, leaving N and acyl C")
    acyl_neighbor = next(iter(acyl_neighbors))
    if topology.atoms[acyl_neighbor].atomic_number != 6:
        raise ValueError("Acyl substituent must be carbon")
    if len(graph[leaving] - {c}) != 2:
        raise ValueError("Expected leaving amide N with two retained substituents")
    graph[c].remove(leaving)
    graph[leaving].remove(c)
    acyl, detached = _component(graph, c), _component(graph, leaving)
    if acyl & detached or (acyl | detached) & {n, og, h}:
        raise ValueError(
            "C-N graph cut must separate substrate from each other and Thr"
        )
    result = source.copy()
    # The nearest intersection of two spheres preserves both the acyl tail and
    # its C-C bond length. Moving only C/O avoids rotating the whole polymer.
    anchor_axis = _unit(source[acyl_neighbor] - source[og])
    separation = float(np.linalg.norm(source[acyl_neighbor] - source[og]))
    cc_length = float(np.linalg.norm(source[c] - source[acyl_neighbor]))
    along = (ester_distance**2 - cc_length**2 + separation**2) / (2 * separation)
    radius_squared = ester_distance**2 - along**2
    if radius_squared <= 1e-8:
        raise ValueError(
            "Fixed OG1/acyl geometry cannot support requested ester distance"
        )
    circle_center = source[og] + along * anchor_axis
    radial = source[c] - circle_center
    radial -= np.dot(radial, anchor_axis) * anchor_axis
    carbon_target = circle_center + np.sqrt(radius_squared) * _unit(radial)
    axis = _unit(source[og] - carbon_target)
    acyl_direction = _unit(source[acyl_neighbor] - carbon_target)
    acyl_angle = float(
        np.degrees(np.arccos(np.clip(np.dot(axis, acyl_direction), -1, 1)))
    )
    if not 95 <= acyl_angle <= 140:
        raise ValueError("Fixed acyl tail requires implausible carbonyl angle")
    acyl_indices = sorted(acyl)
    result[c] = carbon_target
    result[o] = carbon_target - 1.23 * _unit(axis + acyl_direction)
    detached_indices = sorted(detached)
    direction = _unit(source[leaving] - source[c])
    normal = _unit(np.cross(axis, acyl_direction))
    angle = np.radians(leaving_bend_degrees)
    direction = (
        direction * np.cos(angle)
        + np.cross(normal, direction) * np.sin(angle)
        + normal * np.dot(normal, direction) * (1 - np.cos(angle))
    )
    leaving_target = carbon_target + leaving_distance * direction
    result[detached_indices] += leaving_target - source[leaving]
    substituents = sorted(graph[leaving])
    open_valence = -sum(_unit(result[i] - result[leaving]) for i in substituents)
    result[h] = result[leaving] + 1.03 * _unit(open_valence)
    moved = np.flatnonzero(np.linalg.norm(result - source, axis=1) > 1e-10)
    product_graph = [neighbors.copy() for neighbors in graph]
    product_graph[og].remove(h)
    product_graph[h].remove(og)
    for a, b in ((og, c), (leaving, h)):
        product_graph[a].add(b)
        product_graph[b].add(a)
    radii = {1: 1.2, 6: 1.7, 7: 1.55, 8: 1.52, 16: 1.8}
    atom_radii = np.array([radii.get(a.atomic_number, 1.8) for a in topology.atoms])
    orientation = {"axis": None, "degrees": 0.0}
    if orient_leaving:
        result, orientation = _orient_leaving(
            result, sorted(detached | {h}), leaving, product_graph, atom_radii
        )
        moved = np.flatnonzero(np.linalg.norm(result - source, axis=1) > 1e-10)
    contacts = []
    moved_set = set(moved.tolist())
    for i in moved:
        distances = np.linalg.norm(result - result[i], axis=1)
        close = np.flatnonzero(distances < 0.65 * (atom_radii + atom_radii[i]))
        excluded = product_graph[i] | {i}
        for j in close:
            if j in excluded or (j in moved_set and j < i):
                continue
            if any(j in product_graph[k] for k in product_graph[i]):
                continue
            contacts.append(
                {
                    "atoms_1based": [int(i + 1), int(j + 1)],
                    "distance_A": float(distances[j]),
                    "source_distance_A": float(np.linalg.norm(source[i] - source[j])),
                }
            )
    distance_pairs = {
        "OG1_C": (og, c),
        "C_O": (c, o),
        "C_Nleaving": (c, leaving),
        "Nleaving_HG1": (leaving, h),
        "ThrN_H1": (n, ids["thr_h1"]),
        "ThrN_H2": (n, ids["thr_h2"]),
        "ThrN_HG1": (n, h),
    }
    report = {
        "status": "coordinate hypothesis only; not optimized or classically parameterized",
        "atom_count": len(source),
        "source_proton_state": source_proton_state,
        "mapping_1based": vars(mapping).copy(),
        "OG1_C_acyl_angle_degrees": acyl_angle,
        "leaving_rigid_rotation": orientation,
        "leaving_bend_degrees": leaving_bend_degrees,
        "moved_indices_1based": (moved + 1).tolist(),
        "acyl_fragment_indices_1based": [i + 1 for i in acyl_indices],
        "leaving_fragment_indices_1based": [i + 1 for i in detached_indices],
        "distances_A": {
            name: float(np.linalg.norm(result[a] - result[b]))
            for name, (a, b) in distance_pairs.items()
        },
        "close_contacts": contacts,
        "contact_rule": "distance < 0.65 sum VdW radii; product 1-2/1-3 pairs excluded; no PBC",
    }
    return result, report
