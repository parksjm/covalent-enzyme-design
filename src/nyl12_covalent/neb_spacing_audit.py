"""Structured proton diagnostics for source-labeled spacing-test coordinates."""

from __future__ import annotations

from typing import Any

import numpy as np

from .relay_case import RelayCase


def spacing_image_audit(
    model: RelayCase,
    xyz: np.ndarray,
    source_reference: np.ndarray,
) -> dict[str, Any]:
    """Report strict H gates and explicit shared-proton exceptions without hiding them.

    Args:
        model: Unchanged102QM54H explicit case.
        xyz: Full candidate coordinates in source atom order, angstrom.
        source_reference: Frozen native source image1 used for spectator H parents.

    Returns:
        Nearest two heavy atoms for every offending H, spectator-parent changes,
        acid geometry, distances and an explicit strict_H_gate_passed flag.
        Transferring-H-only exceptions are diagnostic, never auto-accepted.
    """
    h = model.hydrogens
    heavy = np.array(
        [i for i in model.qm if model.topology.atoms[i].atomic_number != 1]
    )
    distances = np.linalg.norm(xyz[h, None] - xyz[heavy][None], axis=2)
    reference = np.linalg.norm(
        source_reference[h, None] - source_reference[heavy][None], axis=2
    )
    closest = np.argsort(distances, axis=1)[:, :2]
    nearest = distances[np.arange(len(h)), closest[:, 0]]
    hh = np.linalg.norm(xyz[h, None] - xyz[h][None], axis=2)
    np.fill_diagonal(hh, np.inf)
    m = model.atom_map
    transfer = m["water_h1"] - 1
    changed = []
    offenders = []
    for row, atom in enumerate(h):
        parent = int(heavy[closest[row, 0]])
        if atom != transfer and parent != int(heavy[np.argmin(reference[row])]):
            changed.append(int(atom + 1))
        if nearest[row] < 0.65 or nearest[row] > 1.45:
            offenders.append(
                dict(
                    hydrogen_1based=int(atom + 1),
                    nearest_heavy_1based=parent + 1,
                    nearest_distance_A=float(nearest[row]),
                    second_heavy_1based=int(heavy[closest[row, 1]] + 1),
                    second_distance_A=float(distances[row, closest[row, 1]]),
                )
            )
    acid = float(np.linalg.norm(xyz[m["ash_proton"] - 1] - xyz[m["ash_oxygen"] - 1]))
    no_overlap = bool(nearest.min() >= 0.65 and hh.min() >= 0.65)
    strict = bool(not offenders and not changed and no_overlap and 0.8 <= acid <= 1.3)
    shared_only = bool(
        offenders
        and all(
            row["hydrogen_1based"] == m["water_h1"]
            and {row["nearest_heavy_1based"], row["second_heavy_1based"]}
            == {m["thr_og"], m["leaving_n"]}
            for row in offenders
        )
        and not changed
        and no_overlap
        and 0.8 <= acid <= 1.3
    )
    return dict(
        real_H_count=len(h),
        nearest_H_heavy_min_A=float(nearest.min()),
        nearest_H_heavy_max_A=float(nearest.max()),
        min_H_H_A=float(hh.min()),
        ash_HD2_OD2_A=acid,
        spectator_parent_changes=changed,
        offending_hydrogens=offenders,
        no_compressed_H_contacts=no_overlap,
        strict_H_gate_passed=strict,
        shared_transfer20318_only_requires_root_decision=shared_only,
    )
