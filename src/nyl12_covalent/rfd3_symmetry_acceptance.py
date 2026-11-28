"""Admission checks for independently qualified symmetric RFD3 parents."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

COUPLED_REPEATS = [["A", "E"], ["B", "F"], ["C", "G"], ["D", "H"]]


def require_backbone_symmetry(parent: dict[str, Any]) -> None:
    """Require hash-bound symmetry evidence before admitting an RFD3 parent.

    This consumes a separate geometric qualification; it does not generate
    symmetry operators or certify an unmeasured native output. Missing evidence
    rejects legacy asymmetric parents without rewriting their old receipts.

    Args:
        parent: Canonical RFD3 parent with paths, digests and acceptance receipt.

    Raises:
        ValueError: Evidence is absent, inconsistent or outside declared limits.
    """
    field = "backbone_symmetry_acceptance"
    if field not in parent or field + "_sha256" not in parent:
        raise ValueError("RFD3 backbone symmetry qualification is required")
    path = Path(parent[field])
    if hashlib.sha256(path.read_bytes()).hexdigest() != parent[field + "_sha256"]:
        raise ValueError("RFD3 backbone symmetry receipt hash mismatch")
    receipt = json.loads(path.read_text())
    expected = {
        "schema_version": 1,
        "accepted": True,
        "point_group": "D2",
        "sampler_kind": "symmetry",
        "coupled_repeats": COUPLED_REPEATS,
        "equivalent_mobile_masks": True,
        "operator_closure_accepted": True,
        "fixed_ts_context_retained": True,
        "pdb_path_sha256": parent["pdb_path_sha256"],
        "canonical_mapping_sha256": parent["canonical_mapping_sha256"],
    }
    if any(
        type(receipt.get(key)) is not type(value) or receipt[key] != value
        for key, value in expected.items()
    ):
        raise ValueError("RFD3 backbone symmetry contract is not accepted")
    for metric in ("backbone_rmsd_A", "backbone_max_deviation_A"):
        value = receipt.get(metric)
        limit = receipt.get(metric + "_limit")
        if (
            type(value) not in (float, int)
            or type(limit) not in (float, int)
            or not math.isfinite(value)
            or not math.isfinite(limit)
            or value < 0
            or limit <= 0
            or value > limit
        ):
            raise ValueError(f"RFD3 symmetry residual fails declared limit: {metric}")
