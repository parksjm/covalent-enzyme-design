"""Versioned correction for chemically equivalent Val methyl labels."""

from __future__ import annotations

from typing import Any

import numpy as np

from .rfd3_qualification_audit import coordinate_gates as original_coordinate_gates


def coordinate_gates(
    reference: np.ndarray,
    output: np.ndarray,
    atoms: list[dict[str, Any]],
    native_edges: list[list[int]],
    mobile_residues: list[str],
    forming_pair: list[int],
    thresholds: dict[str, float],
) -> dict[str, Any]:
    """Apply unchanged geometric gates with Val CB correctly nonstereogenic.

    Val has two chemically equivalent methyl substituents. Its CG1/CG2
    permutation does not invert a stereocenter. CA and true Ile/Thr CB checks
    remain unchanged; this wrapper never modifies coordinates or thresholds.
    """
    identities = {f"{a['chain']}{a['residue']}": a["residue_name"] for a in atoms}
    result = original_coordinate_gates(
        reference,
        output,
        atoms,
        native_edges,
        mobile_residues,
        forming_pair,
        thresholds,
    )
    removed, retained = [], []
    for center in result["stereochemical_inversions"]:
        residue, name = center.split(":")
        if name == "CB" and identities[residue] == "VAL":
            removed.append(center)
        else:
            retained.append(center)
    result["stereochemical_inversions"] = retained
    result["equivalent_Val_methyl_label_permutations"] = removed
    result["accepted"] = bool(
        result["fixed_max_A"] <= thresholds["fixed_max_A"]
        and thresholds["mobile_CA_rmsd_min_A"]
        <= result["mobile_CA_rmsd_A"]
        <= thresholds["mobile_CA_rmsd_max_A"]
        and result["peptide_closure_pass"]
        and not result["severe_heavy_clashes"]
        and not retained
    )
    return result
