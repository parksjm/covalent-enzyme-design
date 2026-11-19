"""Bounded deterministic rigid leaving-fragment coordinate search; no engine."""

from __future__ import annotations

from typing import Any

import numpy as np
from parmed import Structure
from scipy.spatial import cKDTree

from .acylenzyme_seed import _component, _unit
from .relay_acyl_seed import RelayAcylAtomMap, _source_graph, build_balanced_acyl_seed


def _rotation(axis: np.ndarray, degrees: float) -> np.ndarray:
    axis = _unit(axis)
    skew = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    angle = np.radians(degrees)
    return np.eye(3) + np.sin(angle) * skew + (1 - np.cos(angle)) * skew @ skew


class _Contacts:
    def __init__(
        self,
        topology: Structure,
        source: np.ndarray,
        moving: np.ndarray,
        graph: list[set[int]],
    ) -> None:
        self.source, self.moving = source, moving
        self.fixed = np.array(sorted(set(range(len(source))) - set(moving.tolist())))
        self.tree = cKDTree(source[self.fixed])
        self.h = np.array([a.atomic_number == 1 for a in topology.atoms])
        self.radii = np.array(
            [
                {1: 1.2, 6: 1.7, 7: 1.55, 8: 1.52}.get(a.atomic_number, 1.8)
                for a in topology.atoms
            ]
        )
        self.excluded = []
        for i in moving:
            neighbors = graph[i] | {int(i)}
            neighbors = neighbors | {j for k in graph[i] for j in graph[k]}
            self.excluded.append(neighbors)

    def assess(self, xyz: np.ndarray) -> tuple[tuple[float, ...], list[dict[str, Any]]]:
        pairs = []
        for row, neighbors in enumerate(
            self.tree.query_ball_point(xyz[self.moving], 2.5)
        ):
            i = int(self.moving[row])
            pairs.extend(
                (i, int(self.fixed[j]))
                for j in neighbors
                if self.fixed[j] not in self.excluded[row]
            )
        for row, i in enumerate(self.moving):
            pairs.extend(
                (int(i), int(j))
                for j in self.moving[row + 1 :]
                if j not in self.excluded[row]
            )
        pairs_array = np.array(pairs, dtype=int)
        if not len(pairs):
            return (0.0, 0.0, 0.0), []
        first, second = pairs_array.T
        distances = np.linalg.norm(xyz[first] - xyz[second], axis=1)
        nh = self.h[first].astype(int) + self.h[second].astype(int)
        severe_limits = np.choose(nh, [1.1, 1.0, 0.8])
        severe = distances < severe_limits
        soft_limits = 0.65 * (self.radii[first] + self.radii[second])
        soft = distances < soft_limits
        score = (
            float(severe.sum()),
            float(np.square(np.maximum(1 - distances / severe_limits, 0)).sum()),
            float(np.square(np.maximum(1 - distances / soft_limits, 0)).sum()),
        )
        contacts = [
            {
                "atoms_1based": [int(first[k] + 1), int(second[k] + 1)],
                "distance_A": float(distances[k]),
                "severe": bool(severe[k]),
                "severe_threshold_A": float(severe_limits[k]),
                "source_distance_A": float(
                    np.linalg.norm(self.source[first[k]] - self.source[second[k]])
                ),
            }
            for k in np.flatnonzero(soft)
        ]
        return score, sorted(contacts, key=lambda row: row["distance_A"])


def search_balanced_acyl_seed(
    topology: Structure,
    source: np.ndarray,
    mapping: RelayAcylAtomMap,
    *,
    active_atom_ids: tuple[int, ...],
    qm_charge: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Search 285 declared rigid leaving-component placements at fixed acyl pose.

    Args:
        topology: Unchanged E-S topology with real source atom identities.
        source: Full finite source coordinates in angstrom.
        mapping: Explicit reactive/water/acid atom identities.
        active_atom_ids: Original one-based active QM/MM atoms allowed to move.
        qm_charge: Unchanged explicit formal QM charge.

    Returns:
        Best coordinate hypothesis and transparent candidate/contact evidence.
        The caller must inspect severe_contact_count; zero does not imply a
        minimum, valid endpoint, or scientific acceptance. Protein and water O
        remain fixed, including original water-H replacement coordinates.

    Raises:
        ValueError: If construction or active-selection prerequisites fail.
    """
    base, construction = build_balanced_acyl_seed(
        topology,
        source,
        mapping,
        active_atom_ids=active_atom_ids,
        qm_charge=qm_charge,
        max_heavy_displacement_A=20.0,
        max_hydrogen_displacement_A=20.0,
        idealize_ester_angle=False,
    )
    graph, ids, _ = _source_graph(topology, source, mapping)
    c, n = ids["carbon"], ids["leaving_n"]
    graph[c].remove(n)
    graph[n].remove(c)
    leaving = _component(graph, n)
    donor = construction["donor_water_h_1based"] - 1
    moving_leave = np.array(sorted(leaving | {donor}), dtype=int)
    moving = np.array(
        sorted(
            set(i - 1 for i in construction["moved_indices_1based"]) | leaving | {donor}
        ),
        dtype=int,
    )
    if not set((moving + 1).tolist()) <= set(active_atom_ids):
        raise ValueError("Search would move atoms outside declared active selection")
    # Restore source graph before applying authoritative product bond edits.
    graph[c].add(n)
    graph[n].add(c)
    for operation, pairs in construction["bond_edits_1based"].items():
        for a1, b1 in pairs:
            a, b = a1 - 1, b1 - 1
            if operation == "remove":
                graph[a].remove(b)
                graph[b].remove(a)
            else:
                graph[a].add(b)
                graph[b].add(a)
    contacts = _Contacts(topology, source, moving, graph)
    direction = _unit(base[n] - base[c])
    u = _unit(np.cross(direction, np.eye(3)[int(np.argmin(np.abs(direction)))]))
    v = np.cross(direction, u)
    retained_c = next(i for i in graph[n] if topology.atoms[i].atomic_number == 6)
    spin_axis = _unit(base[retained_c] - base[n])
    offsets = base[moving_leave] - base[n]
    directions = [(0.0, 0.0, direction)]
    for degrees in (10.0, 20.0, 30.0):
        angle = np.radians(degrees)
        for azimuth in range(0, 360, 60):
            radial = np.cos(np.radians(azimuth)) * u + np.sin(np.radians(azimuth)) * v
            directions.append(
                (
                    degrees,
                    float(azimuth),
                    np.cos(angle) * direction + np.sin(angle) * radial,
                )
            )
    best_score, best, best_contacts, best_choice = None, None, None, None
    rows = []
    for distance in (2.6, 2.8, 3.0):
        for cone, azimuth, vector in directions:
            for spin in (-30.0, -15.0, 0.0, 15.0, 30.0):
                xyz = base.copy()
                xyz[moving_leave] = (
                    offsets @ _rotation(spin_axis, spin).T + base[c] + distance * vector
                )
                score, clashes = contacts.assess(xyz)
                row = dict(
                    distance_A=distance,
                    cone_degrees=cone,
                    azimuth_degrees=azimuth,
                    spin_degrees=spin,
                    score=list(score),
                )
                rows.append(row)
                if best_score is None or score < best_score:
                    best_score, best, best_contacts, best_choice = (
                        score,
                        xyz,
                        clashes,
                        row,
                    )
    if best is None:
        raise ValueError("Declared candidate grid was empty")
    displacement = np.linalg.norm(best - source, axis=1)
    heavy = np.array([a.atomic_number != 1 for a in topology.atoms])
    changed = np.flatnonzero(displacement > 1e-10)
    if not set((changed + 1).tolist()) <= set(active_atom_ids):
        raise ValueError("Selected coordinate moves exceed declared active selection")
    return best, dict(
        status="coordinate search only; no engine or scientific acceptance",
        candidates=rows,
        candidate_count=len(rows),
        best_choice=best_choice,
        severe_contact_count=int(best_score[0]),
        close_contacts=best_contacts,
        moved_indices_1based=(changed + 1).tolist(),
        max_heavy_displacement_A=float(displacement[heavy].max(initial=0)),
        max_hydrogen_displacement_A=float(displacement[~heavy].max(initial=0)),
        scoring="lexicographic severe count, squared severe overlap, squared 0.65-vdW overlap; product 1-2/1-3 excluded; non-PBC",
        severe_thresholds_A=dict(heavy_heavy=1.1, heavy_H=1.0, H_H=0.8),
        water_oxygen_unchanged=bool(
            np.array_equal(best[ids["water_o"]], source[ids["water_o"]])
        ),
        construction=construction,
    )
