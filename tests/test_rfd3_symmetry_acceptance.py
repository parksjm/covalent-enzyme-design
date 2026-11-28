"""Prospective RFD3 admission rejects unqualified or asymmetric parents."""

import json
from pathlib import Path

import pytest

from nyl12_covalent import relay_broad_runtime as runtime
from nyl12_covalent.rfd3_symmetry_acceptance import require_backbone_symmetry

from test_relay_broad_runtime import make_parent


def test_qualified_receipt_and_legacy_hold(tmp_path: Path) -> None:
    """Accept complete synthetic evidence; reject old parents lacking it."""
    parent = make_parent(tmp_path)
    require_backbone_symmetry(parent)
    del parent["backbone_symmetry_acceptance"]
    with pytest.raises(ValueError, match="qualification is required"):
        runtime.verify_parent(parent)


def test_receipt_hash_is_verified(tmp_path: Path) -> None:
    """Changing a receipt invalidates its original digest."""
    parent = make_parent(tmp_path)
    Path(parent["backbone_symmetry_acceptance"]).write_text("{}")
    with pytest.raises(ValueError, match="hash mismatch"):
        require_backbone_symmetry(parent)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("accepted", False),
        ("accepted", 1),
        ("point_group", "C4"),
        ("sampler_kind", "default"),
        ("coupled_repeats", [["A", "B", "C", "D"]]),
        ("equivalent_mobile_masks", False),
        ("operator_closure_accepted", False),
        ("fixed_ts_context_retained", False),
        ("pdb_path_sha256", "other output"),
        ("canonical_mapping_sha256", "other mapping"),
        ("backbone_rmsd_A", None),
        ("backbone_rmsd_A", True),
        ("backbone_rmsd_A", -1),
        ("backbone_rmsd_A", float("nan")),
        ("backbone_rmsd_A", 0.1),
        ("backbone_rmsd_A_limit", 0),
        ("backbone_rmsd_A_limit", float("inf")),
        ("backbone_max_deviation_A", 0.1),
    ],
)
def test_inconsistent_symmetry_is_rejected(
    tmp_path: Path, field: str, value: object
) -> None:
    """Tie approval to the correct output and finite passing measurements."""
    parent = make_parent(tmp_path)
    path = Path(parent["backbone_symmetry_acceptance"])
    receipt = json.loads(path.read_text())
    receipt[field] = value
    path.write_text(json.dumps(receipt))
    parent["backbone_symmetry_acceptance_sha256"] = runtime.digest(path)
    with pytest.raises(ValueError, match="contract|residual"):
        require_backbone_symmetry(parent)
