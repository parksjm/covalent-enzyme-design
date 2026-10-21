"""Bounded geometric feasibility bound for relay-water oxygen sites."""

from __future__ import annotations

from typing import Any

import numpy as np

from .relay_water import _contact_limits

GRID_SETTINGS = {
    "method": "sphere-intersection oxygen grid",
    "n_to_oxygen_A": {"start": 2.6, "stop": 3.4, "step": 0.1, "count": 9},
    "og1_to_oxygen_A": {"start": 2.6, "stop": 3.4, "step": 0.1, "count": 9},
    "azimuth_count": 72,
    "azimuth_endpoint": False,
    "environment": "all topology atoms except WAT1314; no distance cutoff",
}


def rank_oxygen_sites(
    nitrogen: np.ndarray,
    og1: np.ndarray,
    environment: np.ndarray,
    oxygen_limits: np.ndarray,
    *,
    tolerance: float = 0.10,
    retain: int = 10,
    distinct_A: float = 0.4,
) -> dict[str, Any]:
    """Rank oxygen sites on the existing 9x9x72 sphere-intersection grid.

    Args:
        nitrogen: Reactive beta-N coordinate in Angstrom.
        og1: Reactive Thr OG1 coordinate in Angstrom.
        environment: Coordinates of all non-water environment atoms.
        oxygen_limits: Existing per-element oxygen contact limits.
        tolerance: Maximum allowed positive contact shortfall in Angstrom.
        retain: Number of spatially distinct ranked sites to return.
        distinct_A: Minimum separation between retained oxygen sites.

    Returns:
        Eligible-site count and spatially distinct best sites with their worst
        contact shortfalls. This is a finite-grid necessary-condition screen.

    Raises:
        ValueError: If inputs are malformed or the relay axis is degenerate.
    """
    nitrogen = np.asarray(nitrogen, dtype=float)
    og1 = np.asarray(og1, dtype=float)
    environment = np.asarray(environment, dtype=float)
    oxygen_limits = np.asarray(oxygen_limits, dtype=float)
    if nitrogen.shape != (3,) or og1.shape != (3,):
        raise ValueError("N and OG1 coordinates must each have shape (3,)")
    if environment.ndim != 2 or environment.shape[1:] != (3,):
        raise ValueError("Environment coordinates must have shape (n, 3)")
    if oxygen_limits.shape != (len(environment),):
        raise ValueError("Oxygen limits must match the environment atom count")
    if (
        not np.isfinite(nitrogen).all()
        or not np.isfinite(og1).all()
        or not np.isfinite(environment).all()
    ):
        raise ValueError("Coordinates must be finite")
    if not np.isfinite(oxygen_limits).all() or tolerance < 0 or retain < 1:
        raise ValueError("Limits, tolerance, or retain count are invalid")
    axis_vector = og1 - nitrogen
    separation = float(np.linalg.norm(axis_vector))
    if separation < 1.0e-10:
        raise ValueError("Cannot construct an axis from coincident atoms")
    axis = axis_vector / separation
    reference = np.eye(3)[int(np.argmin(np.abs(axis)))]
    first = np.cross(axis, reference)
    first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    sites: list[tuple[float, np.ndarray, float, float]] = []
    for n_distance in np.linspace(2.6, 3.4, 9):
        for og_distance in np.linspace(2.6, 3.4, 9):
            along = (n_distance**2 - og_distance**2 + separation**2) / (2 * separation)
            radius_squared = n_distance**2 - along**2
            if radius_squared <= 0:
                continue
            center = nitrogen + along * axis
            radius = float(np.sqrt(radius_squared))
            for angle in np.linspace(0, 2 * np.pi, 72, endpoint=False):
                point = center + radius * (
                    np.cos(angle) * first + np.sin(angle) * second
                )
                shortfalls = np.maximum(
                    oxygen_limits - np.linalg.norm(environment - point, axis=1), 0.0
                )
                maximum = float(shortfalls.max(initial=0.0))
                sites.append((maximum, point, float(n_distance), float(og_distance)))
    sites.sort(key=lambda item: item[0])
    eligible = sum(site[0] <= tolerance for site in sites)
    distinct: list[tuple[float, np.ndarray, float, float]] = []
    for site in sites:
        if all(np.linalg.norm(site[1] - prior[1]) >= distinct_A for prior in distinct):
            distinct.append(site)
        if len(distinct) == retain:
            break
    return {
        "grid_site_count": len(sites),
        "eligible_site_count": eligible,
        "best_max_shortfall_A": sites[0][0] if sites else None,
        "sites": [
            {
                "rank": rank,
                "oxygen_xyz_A": point.tolist(),
                "n_to_oxygen_A": n_distance,
                "og1_to_oxygen_A": og_distance,
                "maximum_positive_shortfall_A": shortfall,
                "eligible": shortfall <= tolerance,
            }
            for rank, (shortfall, point, n_distance, og_distance) in enumerate(
                distinct, 1
            )
        ],
    }


def _atom_record(atom: Any) -> dict[str, Any]:
    """Return stable atom identity fields for report blockers."""
    residue = atom.residue
    return {
        "atom_index_zero_based": int(atom.idx),
        "residue_index_zero_based": int(residue.idx),
        "residue_index_one_based": int(residue.idx) + 1,
        "residue_name": str(residue.name),
        "atom_name": str(atom.name),
        "element_atomic_number": int(atom.atomic_number),
    }


def _validate_inventory(
    topology: Any, coordinates: np.ndarray
) -> tuple[dict[str, Any], list[int], list[int]]:
    """Validate required residue and proton inventories and identify removals."""
    if len(topology.residues) != 1314:
        raise ValueError("Expected exactly 1314 residues")
    if (
        coordinates.shape != (len(topology.atoms), 3)
        or not np.isfinite(coordinates).all()
    ):
        raise ValueError("Coordinates must be finite and match topology atom count")
    asn, thr, water = (topology.residues[i] for i in (185, 928, 1313))
    if (asn.name, thr.name, water.name) != ("ASN", "THR", "WAT"):
        raise ValueError(
            "Expected ASN186, Thr929, and WAT1314 at validated residue indices"
        )
    asn_atoms = {atom.name: atom for atom in asn.atoms}
    thr_atoms = {atom.name: atom for atom in thr.atoms}
    water_atoms = {atom.name: atom for atom in water.atoms}
    if len(asn_atoms) != len(asn.atoms) or not {"CA", "CB"}.issubset(asn_atoms):
        raise ValueError("ASN186 requires unique CA and CB atoms")
    if asn_atoms["CB"] not in asn_atoms["CA"].bond_partners:
        raise ValueError("ASN186 CB must be bonded to CA")
    if not {"N", "CA", "OG1", "H1", "H2", "H3"}.issubset(thr_atoms):
        raise ValueError("Thr929 must contain N, CA, OG1, and H1/H2/H3")
    if set(water_atoms) != {"O", "H1", "H2"}:
        raise ValueError("WAT1314 must contain O/H1/H2")
    beta_h = [a for a in thr_atoms["N"].bond_partners if a.atomic_number == 1]
    if {a.name for a in beta_h} != {"H1", "H2", "H3"}:
        raise ValueError("Thr929 N must have its three explicitly connected hydrogens")
    sidechain: set[int] = set()
    frontier = [asn_atoms["CB"]]
    while frontier:
        atom = frontier.pop()
        if atom.idx in sidechain:
            continue
        sidechain.add(atom.idx)
        for partner in atom.bond_partners:
            if (
                partner.residue.idx == 185
                and partner.idx != asn_atoms["CA"].idx
                and partner.idx not in sidechain
            ):
                frontier.append(partner)
    sidechain.discard(asn_atoms["CA"].idx)
    removals = sorted(sidechain | {atom.idx for atom in beta_h})
    inventory = {
        "residue_count": len(topology.residues),
        "asn186": {"residue_index_zero_based": 185, "atom_count": len(asn.atoms)},
        "thr929": {
            "residue_index_zero_based": 928,
            "atom_count": len(thr.atoms),
            "beta_n_hydrogens": [a.idx for a in beta_h],
        },
        "wat1314": {
            "residue_index_zero_based": 1313,
            "atom_indices": [water_atoms[n].idx for n in ("O", "H1", "H2")],
        },
        "finite_coordinates": True,
    }
    return inventory, removals, sorted(a.idx for a in water.atoms)


def assess_relay_feasibility(topology: Any, coordinates: np.ndarray) -> dict[str, Any]:
    """Compare original and optimistic fixed-grid oxygen steric feasibility.

    Args:
        topology: Validated ParmEd-like topology for either TI1 target.
        coordinates: Finite coordinates in topology order, Angstrom.

    Returns:
        JSON-serializable inventory, two scenario rankings, blocker identities,
        and interpretation limits. Coordinates and topology are never mutated.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    inventory, removed, water_indices = _validate_inventory(topology, coordinates)
    reactive = {a.name: a.idx for a in topology.residues[928].atoms}
    all_environment = [a.idx for a in topology.atoms if a.idx not in water_indices]
    scenarios = {
        "original_all_atoms": all_environment,
        "optimistic_asn_sidechain_and_beta_n_h_removed": [
            idx for idx in all_environment if idx not in set(removed)
        ],
    }
    results: dict[str, Any] = {}
    for name, indices in scenarios.items():
        limits, _ = _contact_limits(topology, indices)
        ranked = rank_oxygen_sites(
            coordinates[reactive["N"]],
            coordinates[reactive["OG1"]],
            coordinates[indices],
            limits,
        )
        ranked["accepted"] = ranked["eligible_site_count"] > 0
        ranked["acceptance_rule"] = (
            "at least one grid site has maximum shortfall <= 0.10 A"
        )
        for site in ranked["sites"]:
            distances = np.linalg.norm(
                coordinates[indices] - site["oxygen_xyz_A"], axis=1
            )
            shortfalls = np.maximum(limits - distances, 0.0)
            blockers = np.flatnonzero(shortfalls > 0)
            blockers = sorted(blockers, key=lambda pos: shortfalls[pos], reverse=True)
            site["blockers"] = [
                {
                    **_atom_record(topology.atoms[indices[pos]]),
                    "distance_A": float(distances[pos]),
                    "threshold_A": float(limits[pos]),
                    "positive_shortfall_A": float(shortfalls[pos]),
                }
                for pos in blockers
            ]
        results[name] = ranked
    ignored = [_atom_record(topology.atoms[idx]) for idx in removed]
    return {
        "inventory": inventory,
        "grid_settings": GRID_SETTINGS,
        "maximum_allowed_contact_shortfall_A": 0.10,
        "scenarios": results,
        "optimistic_ignored_atoms": ignored,
        "limitations": (
            "A pass is necessary but not sufficient: only oxygen contacts on the finite "
            "existing 9x9x72 grid are tested. It does not move atoms, assess water "
            "hydrogens, side-chain strain, other displaced-atom clashes, energies, "
            "or prove a viable water orientation. A failure bounds this grid only."
        ),
    }
