"""Construct balanced relay acylation guesses; never imply optimized endpoints."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from typing import Any

import numpy as np
from parmed import Structure

from .acylenzyme_seed import AcylAtomMap, _component, _unit, _validated_graph


@dataclass(frozen=True)
class RelayAcylAtomMap:
    """Explicit one-based E-S identities, including unchanged acid states."""

    acyl: AcylAtomMap
    water_o: int
    water_h1: int
    water_h2: int
    ash_oxygen: int
    ash_proton: int
    asp_oxygens: tuple[int, int]


def _source_graph(
    topology: Structure, source: np.ndarray, mapping: RelayAcylAtomMap
) -> tuple[list[set[int]], dict[str, int], list[list[int]]]:
    graph, ids = _validated_graph(topology, mapping.acyl)
    extra = {
        name: getattr(mapping, name) - 1
        for name in ("water_o", "water_h1", "water_h2", "ash_oxygen", "ash_proton")
    }
    asp = [index - 1 for index in mapping.asp_oxygens]
    all_ids = [*ids.values(), *extra.values(), *asp]
    if len(set(all_ids)) != len(all_ids) or any(
        index < 0 or index >= len(source) for index in all_ids
    ):
        raise ValueError("Relay mapping must contain distinct valid atom indices")
    ids.update(extra)
    for key, element in (
        ("water_o", 8),
        ("water_h1", 1),
        ("water_h2", 1),
        ("ash_oxygen", 8),
        ("ash_proton", 1),
    ):
        if topology.atoms[ids[key]].atomic_number != element:
            raise ValueError("Relay atom element mismatch")
    water = {ids["water_h1"], ids["water_h2"]}
    if graph[ids["water_o"]] != water:
        raise ValueError("Mapped source water must contain exactly two hydrogens")
    omitted = []
    h1, h2, ow = ids["water_h1"], ids["water_h2"], ids["water_o"]
    if h2 in graph[h1] or h1 in graph[h2]:
        atoms = [topology.atoms[i] for i in (ow, h1, h2)]
        if (
            any(a.residue is not atoms[0].residue for a in atoms)
            or atoms[0].residue.name != "WAT"
            or [a.name for a in atoms] != ["O", "H1", "H2"]
            or graph[h1] != {ow, h2}
            or graph[h2] != {ow, h1}
        ):
            raise ValueError(
                "Unrecognized water H-H edge; only same-WAT H1-H2 constraint allowed"
            )
        graph[h1].remove(h2)
        graph[h2].remove(h1)
        omitted = [[h1 + 1, h2 + 1]]
    if any(topology.atoms[index].atomic_number != 8 for index in asp):
        raise ValueError("ASP oxygen mapping has incorrect elements")
    required = [
        (ids["thr_og"], ids["proton"]),
        (ids["thr_n"], ids["thr_h1"]),
        (ids["thr_n"], ids["thr_h2"]),
        (ids["water_o"], ids["water_h1"]),
        (ids["water_o"], ids["water_h2"]),
        (ids["ash_oxygen"], ids["ash_proton"]),
    ]
    for parent, hydrogen in required:
        if graph[hydrogen] != {parent} or not (
            0.75 <= np.linalg.norm(source[parent] - source[hydrogen]) <= 1.3
        ):
            raise ValueError("Source proton inventory or parent geometry is invalid")
    for key, expected in (("thr_n", 2), ("ash_oxygen", 1)):
        hydrogens = [i for i in graph[ids[key]] if topology.atoms[i].atomic_number == 1]
        if len(hydrogens) != expected:
            raise ValueError("Source Nalpha/ASH proton inventory is invalid")
    if any(topology.atoms[i].atomic_number == 1 for o in asp for i in graph[o]):
        raise ValueError("ASP272 must remain unprotonated")
    return graph, ids, omitted


def build_balanced_acyl_seed(
    topology: Structure,
    xyz: np.ndarray,
    atommap: RelayAcylAtomMap,
    *,
    active_atom_ids: tuple[int, ...],
    max_heavy_displacement_A: float = 3.0,
    max_hydrogen_displacement_A: float = 5.0,
    qm_charge: int,
    idealize_ester_angle: bool = True,
    ester_distance: float = 1.40,
    leaving_distance: float = 3.00,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build a local ester/amine guess with a label-balanced water relay.

    Args:
        topology: Unmodified E-S topology; never modified or product-parameterized.
        xyz: Full source coordinates in angstrom, in topology atom order.
        atommap: Explicit reactive, relay-water and acid-state identities.
        active_atom_ids: One-based active QM/MM atoms permitted to move.
        max_heavy_displacement_A: Positive heavy-atom movement cap.
        max_hydrogen_displacement_A: Positive hydrogen movement cap.
        qm_charge: Unchanged formal QM charge recorded in the audit.
        idealize_ester_angle: Rotate acyl fragment to 120 degrees if explicitly true.
        ester_distance: Target Thr OG1-carbonyl C separation in angstrom.
        leaving_distance: Target broken carbonyl C-leaving N distance in angstrom.

    Returns:
        Independent coordinate array and JSON-compatible construction evidence.
        No endpoint acceptance, energy or optimization is implied. Original HG1
        replaces the geometrically selected donor water H; that H joins leaving N.

    Raises:
        ValueError: If source chemistry, local geometry or movement bounds fail.
    """
    source = np.asarray(xyz, dtype=float)
    if source.shape != (len(topology.atoms), 3) or not np.isfinite(source).all():
        raise ValueError("Expected finite coordinates with shape (atoms, 3)")
    if any(
        not np.isfinite(cap) or cap <= 0
        for cap in (max_heavy_displacement_A, max_hydrogen_displacement_A)
    ):
        raise ValueError("Displacement cap must be finite and positive")
    if not (1.3 <= ester_distance <= 1.6 and 2.3 <= leaving_distance <= 3.1):
        raise ValueError("Product distances outside explicit seed ranges")
    if not isinstance(qm_charge, int):
        raise ValueError("QM charge must be an explicit integer")
    graph, ids, omitted_constraints = _source_graph(topology, source, atommap)
    c, o, n, og = (ids[key] for key in ("carbon", "oxygen", "leaving_n", "thr_og"))
    acyl = graph[c] - {o, n}
    retained = graph[n] - {c}
    tail = [i for i in retained if topology.atoms[i].atomic_number == 6]
    old_h = [i for i in retained if topology.atoms[i].atomic_number == 1]
    if len(acyl) != 1 or graph[o] != {c} or len(tail) != 1 or len(old_h) != 1:
        raise ValueError("Expected simple amide carbonyl and one-H leaving nitrogen")
    anchor = next(iter(acyl))
    if topology.atoms[anchor].atomic_number != 6:
        raise ValueError("Acyl anchor must be carbon")
    cut = [neighbors.copy() for neighbors in graph]
    cut[c].remove(n)
    cut[n].remove(c)
    acyl_fragment, leaving_fragment = _component(cut, c), _component(cut, n)
    if acyl_fragment & leaving_fragment or (acyl_fragment | leaving_fragment) & {
        ids["thr_n"],
        og,
        ids["proton"],
        ids["water_o"],
    }:
        raise ValueError("Scissile cut must separate substrate components from enzyme")
    result = source.copy()
    target_c = source[og] + ester_distance * _unit(source[c] - source[og])
    axis = _unit(source[og] - target_c)
    old_tail = _unit(source[anchor] - source[c])
    radial = old_tail - np.dot(old_tail, axis) * axis
    if np.linalg.norm(radial) < 1e-8:
        radial = np.cross(axis, np.eye(3)[int(np.argmin(np.abs(axis)))])
    target_tail = -0.5 * axis + np.sqrt(0.75) * _unit(radial)
    cross = np.cross(old_tail, target_tail)
    skew = np.array(
        [[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]]
    )
    cosine = float(np.dot(old_tail, target_tail))
    if cosine < -1 + 1e-8:
        perpendicular = _unit(np.cross(old_tail, axis))
        rotation = 2 * np.outer(perpendicular, perpendicular) - np.eye(3)
    else:
        rotation = np.eye(3) + skew + skew @ skew / (1 + cosine)
    if not idealize_ester_angle:
        rotation, target_tail = np.eye(3), old_tail
    acyl_indices, leaving_indices = sorted(acyl_fragment), sorted(leaving_fragment)
    result[acyl_indices] = (source[acyl_indices] - source[c]) @ rotation.T + target_c
    result[o] = target_c - 1.23 * _unit(axis + target_tail)
    target_n = target_c + leaving_distance * _unit(source[n] - source[c])
    result[leaving_indices] += target_n - source[n]
    water_o = ids["water_o"]
    direction = _unit(source[n] - source[water_o])
    waters = [ids["water_h1"], ids["water_h2"]]
    scores = [
        float(np.dot(_unit(source[h] - source[water_o]), direction)) for h in waters
    ]
    donor = waters[int(np.argmax(scores))]
    result[ids["proton"]] = source[donor]
    result[donor] = result[n] - 1.03 * _unit(
        _unit(result[tail[0]] - result[n]) + _unit(result[old_h[0]] - result[n])
    )
    displacement = np.linalg.norm(result - source, axis=1)
    moved = np.flatnonzero(displacement > 1e-10)
    if not set((moved + 1).tolist()) <= set(active_atom_ids):
        raise ValueError(
            "Construction would move atoms outside explicit active selection"
        )
    caps = np.array(
        [
            max_hydrogen_displacement_A
            if atom.atomic_number == 1
            else max_heavy_displacement_A
            for atom in topology.atoms
        ]
    )
    if np.any(displacement > caps):
        raise ValueError("Construction exceeds declared atom displacement cap")
    edits = {
        "remove": [(og, ids["proton"]), (water_o, donor), (c, n)],
        "add": [(og, c), (water_o, ids["proton"]), (n, donor)],
    }
    product_graph = [neighbors.copy() for neighbors in graph]
    for operation, pairs in edits.items():
        for a, b in pairs:
            if operation == "remove":
                product_graph[a].remove(b)
                product_graph[b].remove(a)
            else:
                product_graph[a].add(b)
                product_graph[b].add(a)
    counts = {
        key: sum(topology.atoms[i].atomic_number == 1 for i in product_graph[ids[key]])
        for key in ("thr_n", "thr_og", "leaving_n", "water_o", "ash_oxygen")
    }
    expected = dict(thr_n=2, thr_og=0, leaving_n=2, water_o=2, ash_oxygen=1)
    if counts != expected:
        raise ValueError("Product proton inventory failed")
    contacts = []
    radii = np.array(
        [
            {1: 1.2, 6: 1.7, 7: 1.55, 8: 1.52}.get(a.atomic_number, 1.8)
            for a in topology.atoms
        ]
    )
    for i in moved:
        distances = np.linalg.norm(result - result[i], axis=1)
        excluded = product_graph[i] | {int(i)}
        excluded.update(
            j for neighbor in product_graph[i] for j in product_graph[neighbor]
        )
        for j in np.flatnonzero(distances < 0.65 * (radii[i] + radii)):
            if j in excluded or (j in moved and j < i):
                continue
            contacts.append(
                {
                    "atoms_1based": [int(i + 1), int(j + 1)],
                    "distance_A": float(distances[j]),
                }
            )
    return result, {
        "status": "coordinate hypothesis only; unbiased QM/MM relaxation required",
        "source_coordinate_sha256": hashlib.sha256(
            source.astype("<f8").tobytes()
        ).hexdigest(),
        "mapping_1based": asdict(atommap),
        "omitted_topology_constraints_1based": omitted_constraints,
        "constraint_policy": "same-WAT O/H1/H2 triangle: omit H1-H2 only from chemical graph; topology unchanged",
        "donor_water_h_1based": donor + 1,
        "donor_selection_cosines": scores,
        "idealize_ester_angle": idealize_ester_angle,
        "ester_angle_degrees": float(
            np.degrees(np.arccos(np.clip(np.dot(axis, target_tail), -1, 1)))
        ),
        "geometry_steps": [
            "rigid acyl translation; optional explicit 120 degree rotation",
            "planar carbonyl O at 1.23 A",
            "rigid leaving component translation to explicit C-N target",
            "relay H label permutation; original water O/spectator H fixed",
        ],
        "bond_edits_1based": {
            key: [[a + 1, b + 1] for a, b in pairs] for key, pairs in edits.items()
        },
        "moved_indices_1based": (moved + 1).tolist(),
        "max_displacement_A": float(displacement.max(initial=0)),
        "displacements_A": {str(i + 1): float(displacement[i]) for i in moved},
        "close_contacts": contacts,
        "contact_rule": "non-PBC distance < 0.65 sum vdW radii; omit product 1-2/1-3",
        "proton_inventory": counts,
        "asp_protons": 0,
        "atom_count": len(source),
        "atom_order_and_elements_unchanged": True,
        "hydrogen_count": sum(a.atomic_number == 1 for a in topology.atoms),
        "mass_inventory_unchanged": True,
        "water_oxygen_unchanged": True,
        "qm_charge_before": qm_charge,
        "qm_charge_after": qm_charge,
        "classical_product_topology_valid": False,
        "distances_A": {
            "OG1_C": float(np.linalg.norm(result[og] - result[c])),
            "C_O": float(np.linalg.norm(result[c] - result[o])),
            "C_Nleaving": float(np.linalg.norm(result[c] - result[n])),
        },
    }
