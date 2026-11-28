"""Bounded explicit-hydrogen relay paths; structural hypotheses, not mechanisms."""

from __future__ import annotations

from typing import Any

import numpy as np

from nyl12_covalent.active_site_geometry import projected_herschlag_phi
from nyl12_covalent.nalpha_water_relay import angle, hbond, minimum_image


def best_edge(
    xyz: np.ndarray, donor: dict, acceptor: int, cell: np.ndarray, criteria: dict
) -> dict[str, Any]:
    """Measure the best topology-owned donor hydrogen against one acceptor.

    Args:
        xyz: Frame coordinates in angstroms.
        donor: Resolved chemical role containing hydrogen indices.
        acceptor: Zero-based acceptor index.
        cell: Periodic lattice matrix.
        criteria: da_A, ha_A and angle_deg thresholds.

    Returns:
        Best hydrogen geometry and continuous threshold violation.
    """
    edges = []
    for hydrogen in donor["h_indices"]:
        result = dict(hbond(xyz, donor["index"], hydrogen, acceptor, cell))
        result["pass"] = (
            result["da_A"] <= criteria["da_A"]
            and result["ha_A"] <= criteria["ha_A"]
            and result["angle_deg"] >= criteria["angle_deg"]
        )
        result["violation"] = (
            max(0.0, result["da_A"] - criteria["da_A"])
            + max(0.0, result["ha_A"] - criteria["ha_A"])
            + max(0.0, criteria["angle_deg"] - result["angle_deg"]) / 60
        )
        edges.append({**result, "hydrogen_index1": int(hydrogen + 1)})
    if not edges:
        raise ValueError("Donor has no explicitly bonded hydrogens")
    return min(edges, key=lambda edge: (edge["violation"], edge["ha_A"]))


def screen_paths(
    xyz: np.ndarray, roles: list[dict], target: int, cell: np.ndarray, criteria: dict
) -> list[dict]:
    """Find directed routes with at most two waters and one declared internal role.

    Prefix heavy-atom distances must pass. Explicit hydrogen failures remain
    nearpass diagnostics; strict routes require every D-H-A edge to pass. N is a conditional reaction-prep
    acceptor; this does not assert equilibrium amide basicity.

    Args:
        xyz: Frame coordinates.
        roles: Topology-validated donor and acceptor roles.
        target: Leaving amide nitrogen index.
        cell: Periodic lattice.
        criteria: Hydrogen-bond thresholds.

    Returns:
        Strict routes plus best final-edge nearpass for each starting donor.
    """
    starts = [role for role in roles if role.get("start")]
    # Any strict route has at most four edges (2 waters + 1 internal role).
    radius = 4 * criteria["da_A"]
    internal = [role for role in roles if role.get("internal") and role["acceptor"]]
    distances = np.linalg.norm(
        minimum_image(
            np.array([xyz[r["index"]] for r in internal]).reshape(-1, 3) - xyz[target],
            cell,
        ),
        axis=1,
    )
    local = [
        role
        for role, distance in zip(internal, distances, strict=True)
        if distance <= radius
    ]
    local_xyz = np.array([xyz[r["index"]] for r in local]).reshape(-1, 3)
    neighbors: dict[int, list[dict]] = {}

    def nearby(donor: dict) -> list[dict]:
        index = donor["index"]
        if index not in neighbors:
            distances = np.linalg.norm(
                minimum_image(local_xyz - xyz[index], cell), axis=1
            )
            neighbors[index] = [
                r
                for r, d in zip(local, distances, strict=True)
                if d <= criteria["da_A"]
            ]
        return neighbors[index]

    cache: dict[tuple[int, int], dict] = {}

    def edge(donor: dict, acceptor: int) -> dict:
        key = donor["index"], acceptor
        if key not in cache:
            cache[key] = best_edge(xyz, donor, acceptor, cell, criteria)
        return cache[key]

    results = []
    for start in starts:
        found = []

        def visit(path: list[dict], legs: list[dict], nw: int, ni: int) -> None:
            donor = path[-1]
            last = edge(donor, target)
            full = [*legs, last]
            found.append(
                {
                    "start": start["label"],
                    "path": [r["label"] for r in path] + ["leaving_amide_N"],
                    "atom_indices1": [int(r["index"] + 1) for r in path]
                    + [int(target + 1)],
                    "edges": full,
                    "water_count": nw,
                    "internal_count": ni,
                    "strict": all(bool(e["pass"]) for e in full),
                    "violation": float(sum(e["violation"] for e in full)),
                    "states": [r["state"] for r in path],
                    "conditional_roles": [r["label"] for r in path if r["conditional"]],
                    "leaving_N_acceptor_conditional": True,
                }
            )
            for nxt in nearby(donor):
                if nxt["index"] in {r["index"] for r in path} or not nxt["h_indices"]:
                    continue
                new_nw, new_ni = nw + int(nxt["water"]), ni + int(not nxt["water"])
                if new_nw > 2 or new_ni > 1:
                    continue
                measured = edge(donor, nxt["index"])
                visit([*path, nxt], [*legs, measured], new_nw, new_ni)

        visit([start], [], 0, 0)
        strict = [row for row in found if row["strict"]]
        results.extend(strict or [min(found, key=lambda row: row["violation"])])
    return results


def catalytic_geometry(
    xyz: np.ndarray, mapped: dict[str, int], cell: np.ndarray
) -> dict[str, float]:
    """Measure fresh PBC-aware catalytic geometry with existing signed phi math.

    Args:
        xyz: Frame coordinates.
        mapped: og1, nalpha, carbon, oxygen, nleave, adjacent_c, oxyanion1/2.
        cell: Periodic lattice.

    Returns:
        d1/d2/d3, attack angle and signed projected Herschlag phi.
    """
    points = {
        key: minimum_image(xyz[index] - xyz[mapped["carbon"]], cell)
        for key, index in mapped.items()
    }
    return {
        "d1": float(np.linalg.norm(points["og1"])),
        "d2": float(np.linalg.norm(points["oxygen"] - points["oxyanion1"])),
        "d3": float(np.linalg.norm(points["oxygen"] - points["oxyanion2"])),
        "attack": angle(points["og1"], points["carbon"], points["oxygen"]),
        "projected_phi": projected_herschlag_phi(
            *[
                tuple(points[k])
                for k in ("og1", "oxygen", "carbon", "nleave", "adjacent_c")
            ]
        ),
    }
