"""Select an explicitly configured phi sign without changing geometry tiers."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from itertools import permutations
import math
from typing import Any


def _pose_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        row["tier"] == "B",
        float(row["d1"]),
        str(row["carbonyl_c"]),
        str(row["carbonyl_o"]),
        str(row["leaving_n"]),
    )


def evaluate_signed_amide_gate(
    rows: Sequence[Mapping[str, Any]],
    *,
    target_phi_sign: int,
    sites: Sequence[str] = ("AE", "BF", "CG", "DH"),
    ligands: Sequence[str] = ("I", "J", "K", "L"),
    amides_per_ligand: int = 3,
) -> dict[str, Any]:
    """Filter every amide by an explicit phi sign before Tier A/B assignment.

    Args:
        rows: Complete, chemically audited site/ligand/amide geometry rows.
        target_phi_sign: Required +1 or -1 under the unchanged phi definition.
        sites: Ordered physical catalytic-site labels.
        ligands: Distinct ligand identities, one for each occupied site.
        amides_per_ligand: Number of graph-identified internal amides per ligand.

    Returns:
        All eligible rows, a maximum-cardinality one-to-one assignment in site
        order (unfilled slots are None), counts and the explicit policy.
        Input rows and their legacy tiers are not modified.

    Raises:
        ValueError: Coverage, identity, tier or numerical geometry is invalid.
    """
    if type(target_phi_sign) is not int or target_phi_sign not in (-1, 1):
        raise ValueError("target_phi_sign must be integer +1 or -1; bool is invalid")
    if (
        not sites
        or len(sites) != len(ligands)
        or len(set(sites)) != len(sites)
        or len(set(ligands)) != len(ligands)
        or amides_per_ligand < 1
    ):
        raise ValueError("Expected unique, equally sized site/ligand sets and amides")
    expected = len(sites) * len(ligands) * amides_per_ligand
    if len(rows) != expected:
        raise ValueError(f"Expected all {expected} site/ligand/amide rows")
    seen: set[tuple[Any, ...]] = set()
    coverage: Counter[tuple[str, str]] = Counter()
    identities: dict[str, set[tuple[str, str, str]]] = {lig: set() for lig in ligands}
    eligible = []
    for row in rows:
        site, lig = row["site"], row["ligand"]
        if site not in sites or lig not in ligands:
            raise ValueError("Unexpected physical site or ligand")
        amide = tuple(str(row[k]) for k in ("carbonyl_c", "carbonyl_o", "leaving_n"))
        identity = (site, lig, *amide)
        if any(not name for name in amide) or identity in seen:
            raise ValueError("Empty or duplicate amide identity")
        seen.add(identity)
        coverage[site, lig] += 1
        identities[lig].add(amide)
        phi, d1 = float(row["herschlag_phi"]), float(row["d1"])
        if not math.isfinite(phi) or not math.isfinite(d1) or d1 <= 0:
            raise ValueError("Expected finite phi and positive finite distance")
        if row["tier"] not in {"A", "B", "C", "reject"}:
            raise ValueError("Unknown legacy tier")
        if row["tier"] in {"A", "B"} and phi * target_phi_sign > 0:
            eligible.append(dict(row))
    if any(
        coverage[site, lig] != amides_per_ligand for site in sites for lig in ligands
    ):
        raise ValueError("Incomplete site/ligand coverage")
    if any(len(value) != amides_per_ligand for value in identities.values()):
        raise ValueError("Inconsistent graph-amide identities across sites")
    best: dict[tuple[str, str], dict[str, Any]] = {}
    for row in eligible:
        pair = (row["site"], row["ligand"])
        if pair not in best or _pose_key(row) < _pose_key(best[pair]):
            best[pair] = row
    choices = []
    for permutation in permutations(ligands):
        assignment = [
            best.get((site, lig)) for site, lig in zip(sites, permutation, strict=True)
        ]
        present = [row for row in assignment if row is not None]
        key = (
            -len(present),
            sum(row["tier"] == "B" for row in present),
            sum(float(row["d1"]) for row in present),
            permutation,
        )
        choices.append((key, assignment))
    _, assigned = min(choices, key=lambda value: value[0])
    assigned_count = sum(row is not None for row in assigned)
    return {
        "policy": "explicit_signed_phi_AB_before_max_cardinality_then_tier_total_d1",
        "target_phi_sign": target_phi_sign,
        "input_row_count": len(rows),
        "legacy_tiers_modified": False,
        "eligible_ab": sorted(
            eligible, key=lambda row: (row["site"], row["ligand"], _pose_key(row))
        ),
        "eligible_ab_count": len(eligible),
        "assigned_site_order": list(sites),
        "assigned": assigned,
        "assigned_ab_count": assigned_count,
        "passes_pose_gate": assigned_count > 0,
    }


def evaluate_negative_amide_gate(
    rows: Sequence[Mapping[str, Any]],
    *,
    sites: Sequence[str] = ("AE", "BF", "CG", "DH"),
    ligands: Sequence[str] = ("I", "J", "K", "L"),
    amides_per_ligand: int = 3,
) -> dict[str, Any]:
    """Retain historical negative-sign diagnostics, not campaign acceptance.

    Args:
        rows: Complete chemically audited geometry rows.
        sites: Ordered physical catalytic-site labels.
        ligands: Distinct ligand identities, one per occupied site.
        amides_per_ligand: Number of graph-identified internal amides per ligand.

    Returns:
        The negative-sign diagnostic with historical policy and count aliases.
    """
    result = evaluate_signed_amide_gate(
        rows,
        target_phi_sign=-1,
        sites=sites,
        ligands=ligands,
        amides_per_ligand=amides_per_ligand,
    )
    return {
        **result,
        "policy": "negative_phi_AB_before_max_cardinality_then_tier_total_d1",
        "eligible_negative_ab": result["eligible_ab"],
        "eligible_negative_ab_count": result["eligible_ab_count"],
        "assigned_negative_ab_count": result["assigned_ab_count"],
    }
