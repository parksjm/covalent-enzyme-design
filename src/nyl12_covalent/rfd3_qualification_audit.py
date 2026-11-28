"""Numerical, source-mapped structural gates for fixed-length RFD3 qualification."""

from __future__ import annotations

from typing import Any

import numpy as np


def align_on_fixed(
    reference: np.ndarray, mobile: np.ndarray, fixed: np.ndarray
) -> np.ndarray:
    """Align output coordinates onto immutable source atoms using a proper rotation."""
    if (
        reference.shape != mobile.shape
        or reference.ndim != 2
        or reference.shape[1] != 3
    ):
        raise ValueError("Output coordinate shape differs from source")
    if (
        not np.isfinite(reference).all()
        or not np.isfinite(mobile).all()
        or fixed.sum() < 3
    ):
        raise ValueError(
            "Invalid coordinates or insufficient immutable alignment anchors"
        )
    x, y = mobile[fixed], reference[fixed]
    cx, cy = x.mean(axis=0), y.mean(axis=0)
    u, _, vt = np.linalg.svd((x - cx).T @ (y - cy))
    parity = np.eye(3)
    parity[2, 2] = np.linalg.det(u @ vt)
    return (mobile - cx) @ (u @ parity @ vt) + cy


def coordinate_gates(
    reference: np.ndarray,
    output: np.ndarray,
    atoms: list[dict[str, Any]],
    native_edges: list[list[int]],
    mobile_residues: list[str],
    forming_pair: list[int],
    thresholds: dict[str, float],
) -> dict[str, Any]:
    """Check preservation, local movement, peptide closure, chirality and clashes.

    Args:
        reference: Source heavy coordinates in declared atom order.
        output: Generated coordinates mapped into exactly that order.
        atoms: Source chain, local residue and atom-name records.
        native_edges: Source heavy bond pairs, without an invented forming bond.
        mobile_residues: Explicit occupied-repeat positions permitted to move.
        forming_pair: The mapped TS-like Thr OG1/carbonyl-C contact.
        thresholds: Predeclared numerical geometry limits.

    Returns:
        Acceptance and independent numerical gate records.
    """
    if len(atoms) != len(reference):
        raise ValueError("Atom map and source coordinate counts differ")
    keys = [(a["chain"], int(a["residue"]), a["name"]) for a in atoms]
    if len(set(keys)) != len(keys):
        raise ValueError("Ambiguous source atom identities")
    lookup = {key: index for index, key in enumerate(keys)}
    mobile_set = set(mobile_residues)
    moving = np.array([f"{c}{p}" in mobile_set for c, p, _ in keys])
    if not moving.any():
        raise ValueError("Empty mobile residue mask")
    fixed = ~moving
    aligned = align_on_fixed(reference, output, fixed)
    displacement = np.linalg.norm(aligned - reference, axis=1)
    ca = moving & np.array([name == "CA" for _, _, name in keys])
    if not ca.any():
        raise ValueError("Mobile mask has no backbone CA atoms")
    rmsd = float(np.sqrt(np.mean(displacement[ca] ** 2)))
    closures = []
    stereochemistry = []
    for (chain, position, name), index in lookup.items():
        if name == "C" and (chain, position + 1, "N") in lookup:
            other = lookup[(chain, position + 1, "N")]
            if moving[index] or moving[other]:
                distance = float(np.linalg.norm(aligned[index] - aligned[other]))
                closures.append(
                    {"chain": chain, "residue": position, "distance_A": distance}
                )
        if name == "CA" and (chain, position, "CB") in lookup:
            centers = [("CA", ("N", "C", "CB"))]
            if (chain, position, "OG1") in lookup:
                centers.append(("CB", ("CA", "OG1", "CG2")))
            if (chain, position, "CG1") in lookup and (
                chain,
                position,
                "CG2",
            ) in lookup:
                centers.append(("CB", ("CA", "CG1", "CG2")))
            for center_name, neighbors in centers:
                center = lookup[(chain, position, center_name)]
                indices = [lookup[(chain, position, n)] for n in neighbors]
                before = float(np.linalg.det(reference[indices] - reference[center]))
                after = float(np.linalg.det(aligned[indices] - aligned[center]))
                if before * after <= 0:
                    stereochemistry.append(f"{chain}{position}:{center_name}")
    adjacency: list[set[int]] = [set() for _ in atoms]
    for a, b in native_edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    excluded = {tuple(sorted(forming_pair))}
    for index, neighbors in enumerate(adjacency):
        for other in neighbors | set().union(*(adjacency[n] for n in neighbors)):
            excluded.add(tuple(sorted((index, other))))
    clashes = []
    for index in np.flatnonzero(moving):
        distances = np.linalg.norm(aligned - aligned[index], axis=1)
        for other in np.flatnonzero(distances < thresholds["clash_min_A"]):
            if (
                index != other
                and tuple(sorted((int(index), int(other)))) not in excluded
            ):
                if moving[other] and other < index:
                    continue
                clashes.append(
                    {
                        "atom_indices": [int(index), int(other)],
                        "distance_A": float(distances[other]),
                    }
                )
    fixed_max = float(displacement[fixed].max())
    closure_pass = all(
        thresholds["peptide_min_A"] <= x["distance_A"] <= thresholds["peptide_max_A"]
        for x in closures
    )
    passed = (
        fixed_max <= thresholds["fixed_max_A"]
        and thresholds["mobile_CA_rmsd_min_A"]
        <= rmsd
        <= thresholds["mobile_CA_rmsd_max_A"]
        and closure_pass
        and not clashes
        and not stereochemistry
    )
    return {
        "accepted": bool(passed),
        "fixed_max_A": fixed_max,
        "mobile_CA_rmsd_A": rmsd,
        "peptide_closure_pass": closure_pass,
        "peptide_closures": closures,
        "stereochemical_inversions": stereochemistry,
        "severe_heavy_clashes": clashes,
        "forming_contact_excluded_from_clash_gate": forming_pair,
        "thresholds": thresholds,
    }
