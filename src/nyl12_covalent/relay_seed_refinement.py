"""Refine a coordinate-only acyl guess by one declared internal torsion."""

from __future__ import annotations

from typing import Any

import numpy as np
from parmed import Structure

from .acylenzyme_seed import _component, _unit
from .relay_acyl_seed import RelayAcylAtomMap, _source_graph
from .relay_seed_search import _Contacts, _rotation


def refine_acyl_torsion_seed(
    topology: Structure,
    source: np.ndarray,
    mapping: RelayAcylAtomMap,
    base_xyz: np.ndarray,
    base_evidence: dict[str, Any],
    *,
    active_atom_ids: tuple[int, ...],
    max_heavy_displacement_A: float = 3.0,
    max_hydrogen_displacement_A: float = 5.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Select among 24 acyl-arm rotations, keeping ester/leaving/water fixed.

    Args:
        topology: Unchanged E-S topology defining the real acyl-adjacent bond.
        source: Original full-system coordinates in angstrom.
        mapping: Explicit reactive, water and acid identities.
        base_xyz: Best output of the bounded leaving-fragment coordinate search.
        base_evidence: Full evidence returned by search_balanced_acyl_seed.
        active_atom_ids: Original one-based active QM/MM atom identities.
        max_heavy_displacement_A: Explicit final heavy-atom displacement cap.
        max_hydrogen_displacement_A: Explicit final hydrogen displacement cap.

    Returns:
        Best torsion coordinates and complete construction evidence, including
        all 24 scores. A nonzero severe_contact_count means the caller must
        reject the guess. Zero severe contacts never establishes optimization.

    Raises:
        ValueError: If geometry, source graph, active region or final caps fail.
    """
    source, base = np.asarray(source, dtype=float), np.asarray(base_xyz, dtype=float)
    if (
        source.shape != (len(topology.atoms), 3)
        or base.shape != source.shape
        or not (np.isfinite(source).all() and np.isfinite(base).all())
    ):
        raise ValueError("Expected finite matching full-system source/base coordinates")
    caps = (max_heavy_displacement_A, max_hydrogen_displacement_A)
    if any(not np.isfinite(cap) or cap <= 0 for cap in caps):
        raise ValueError("Final displacement caps must be positive and finite")
    graph, ids, _ = _source_graph(topology, source, mapping)
    c, o, n = ids["carbon"], ids["oxygen"], ids["leaving_n"]
    neighbors = graph[c] - {o, n}
    if len(neighbors) != 1:
        raise ValueError("Expected one actual carbonyl-adjacent acyl atom")
    anchor = next(iter(neighbors))
    graph[c].remove(anchor)
    graph[anchor].remove(c)
    arm = np.array(sorted(_component(graph, anchor)), dtype=int)
    if {c, o, n, ids["thr_og"], ids["water_o"]} & set(arm.tolist()):
        raise ValueError("Acyl torsion cut failed to isolate acyl arm")
    graph[c].add(anchor)
    graph[anchor].add(c)
    construction = dict(base_evidence["construction"])
    for operation, pairs in construction["bond_edits_1based"].items():
        for first, second in pairs:
            a, b = first - 1, second - 1
            if operation == "remove":
                graph[a].remove(b)
                graph[b].remove(a)
            else:
                graph[a].add(b)
                graph[b].add(a)
    initial_moved = np.flatnonzero(np.linalg.norm(base - source, axis=1) > 1e-10)
    moving = np.array(sorted(set(initial_moved.tolist()) | set(arm.tolist())))
    if not set((moving + 1).tolist()) <= set(active_atom_ids):
        raise ValueError("Acyl torsion would move atoms outside original active region")
    scorer = _Contacts(topology, source, moving, graph)
    axis = _unit(base[anchor] - base[c])
    best_score, best, chosen, contacts = None, None, None, None
    rows = []
    for angle in range(-180, 180, 15):
        candidate = base.copy()
        candidate[arm] = (base[arm] - base[c]) @ _rotation(axis, angle).T + base[c]
        candidate[anchor] = base[anchor]
        score, collisions = scorer.assess(candidate)
        rows.append({"torsion_degrees": angle, "score": list(score)})
        if best_score is None or score < best_score:
            best_score, best, chosen, contacts = score, candidate, angle, collisions
    if best is None:
        raise ValueError("Declared torsion grid was empty")
    displacement = np.linalg.norm(best - source, axis=1)
    heavy = np.array([atom.atomic_number != 1 for atom in topology.atoms])
    heavy_max = float(displacement[heavy].max(initial=0))
    hydrogen_max = float(displacement[~heavy].max(initial=0))
    if (
        heavy_max > max_heavy_displacement_A
        or hydrogen_max > max_hydrogen_displacement_A
    ):
        raise ValueError(
            "Selected torsion seed exceeds conservative final movement cap"
        )
    summary = {
        "candidate_count": len(rows),
        "candidates": rows,
        "torsion_axis_1based": [c + 1, anchor + 1],
        "torsion_arm_1based": (arm + 1).tolist(),
        "best_torsion_degrees": chosen,
        "severe_contact_count": int(best_score[0]),
        "previous_leaving_choice": base_evidence["best_choice"],
        "max_heavy_displacement_A": heavy_max,
        "max_hydrogen_displacement_A": hydrogen_max,
        "final_caps_A": {
            "heavy": max_heavy_displacement_A,
            "hydrogen": max_hydrogen_displacement_A,
        },
    }
    construction.update(
        status="coordinate hypothesis only; unbiased QM/MM relaxation required",
        search={
            "leaving": {k: v for k, v in base_evidence.items() if k != "construction"},
            "acyl_torsion": summary,
        },
        severe_contact_count=int(best_score[0]),
        close_contacts=contacts,
        moved_indices_1based=(np.flatnonzero(displacement > 1e-10) + 1).tolist(),
        max_displacement_A=float(displacement.max(initial=0)),
        max_heavy_displacement_A=heavy_max,
        max_hydrogen_displacement_A=hydrogen_max,
        displacements_A={
            str(i + 1): float(displacement[i])
            for i in np.flatnonzero(displacement > 1e-10)
        },
    )
    construction["distances_A"] = {
        key: float(np.linalg.norm(best[a] - best[b]))
        for key, a, b in (
            ("OG1_C", ids["thr_og"], c),
            ("C_O", c, o),
            ("C_Nleaving", c, n),
        )
    }
    return best, construction
