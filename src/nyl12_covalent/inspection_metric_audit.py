"""Separate coordinate serialization precision from reversible conversion checks."""

from __future__ import annotations

from typing import Any

from .catalytic_alignment_geometry import wrapped_delta
from .negative_amide_gate import evaluate_signed_amide_gate

METRICS = ("d1", "d2", "d3", "a_bd", "herschlag_phi")


def serialization_metric_audit(
    native: dict[str, Any], serialized: dict[str, Any], restored: dict[str, Any]
) -> dict[str, Any]:
    """Audit all 48 rows without treating PDB quantization as a mapping change.

    Args:
        native: Full-precision source-coordinate geometry evaluation.
        serialized: Geometry of the actual PDB supplied to the converter.
        restored: Geometry after reversing the converter's atom/residue mapping.

    Returns:
        Exact conversion errors and explicit per-row precision/gate diagnostics.
    """

    def index(evaluation: dict[str, Any]) -> dict[tuple, dict]:
        rows = evaluation["all_amide_site_rows"]
        result = {(r["site"], r["ligand"], r["carbonyl_c"]): r for r in rows}
        if len(rows) != 48 or len(result) != 48:
            raise ValueError("Expected all 48 unique amide/site rows")
        return result

    def error(key: str, first: float, second: float) -> float:
        return (
            abs(wrapped_delta(first, second))
            if key == "herschlag_phi"
            else abs(first - second)
        )

    raw, pdb, back = index(native), index(serialized), index(restored)
    if raw.keys() != pdb.keys() or raw.keys() != back.keys():
        raise ValueError("Serialization metric row identities differ")
    conversion = {k: max(error(k, pdb[i][k], back[i][k]) for i in raw) for k in METRICS}
    if any(
        v > (0.2 if k in ("a_bd", "herschlag_phi") else 0.003)
        for k, v in conversion.items()
    ):
        raise ValueError("Roundtrip metric mismatch")
    assigned_raw = {
        (r["site"], r["ligand"], r["carbonyl_c"]) for r in native["assigned"]
    }
    assigned_pdb = {
        (r["site"], r["ligand"], r["carbonyl_c"]) for r in serialized["assigned"]
    }
    diagnostics = [
        {
            "site": i[0],
            "ligand": i[1],
            "carbonyl_c": i[2],
            "native_metrics": {k: raw[i][k] for k in METRICS},
            "serialized_metrics": {k: pdb[i][k] for k in METRICS},
            "metric_errors": {k: error(k, raw[i][k], pdb[i][k]) for k in METRICS},
            "native_tier": raw[i]["tier"],
            "serialized_tier": pdb[i]["tier"],
            "native_legacy_assigned": i in assigned_raw,
            "serialized_legacy_assigned": i in assigned_pdb,
            "positive_AB_membership_changed": (
                raw[i]["tier"] in ("A", "B") and raw[i]["herschlag_phi"] > 0
            )
            != (pdb[i]["tier"] in ("A", "B") and pdb[i]["herschlag_phi"] > 0),
        }
        for i in raw
    ]
    return {
        "metric_max_errors": conversion,
        "roundtrip_metric_comparison_basis": "serialized PDB supplied to converter versus restored PDB; unchanged thresholds",
        "native_to_PDB_quantization_metric_max_errors": {
            k: max(d["metric_errors"][k] for d in diagnostics) for k in METRICS
        },
        "native_to_PDB_quantization_all48": diagnostics,
        "native_to_PDB_tier_changes": [
            d for d in diagnostics if d["native_tier"] != d["serialized_tier"]
        ],
        "native_to_PDB_positive_AB_membership_changes": [
            d for d in diagnostics if d["positive_AB_membership_changed"]
        ],
        "native_to_PDB_legacy_assignment_changed": assigned_raw != assigned_pdb,
        "native_positive_phi_AB_gate": evaluate_signed_amide_gate(
            native["all_amide_site_rows"], target_phi_sign=1
        ),
        "serialized_positive_phi_AB_gate": evaluate_signed_amide_gate(
            serialized["all_amide_site_rows"], target_phi_sign=1
        ),
    }
