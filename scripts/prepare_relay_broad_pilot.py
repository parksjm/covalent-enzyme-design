#!/usr/bin/env python3
"""Freeze the broad LigandMPNN pilot mask and bounded queue request."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.relay_broad_pilot import (  # noqa: E402
    protection_census,
    sha256_file,
    symmetry_manifest,
)

CAMPAIGN = ROOT / "design/nyl12_relay_ts_ensemble_20261126"
ORIGINAL_PARENTS = {
    "medoid_01_delivery_full_010_frame_068": CAMPAIGN
    / "templates/medoid_01_delivery_full_010_frame_068/complex_heavy.pdb",
    "medoid_03_delivery_full_007_frame_095": CAMPAIGN
    / "templates/medoid_03_delivery_full_007_frame_095/complex_heavy.pdb",
}
CENSUS_PARENTS = [
    CAMPAIGN / f"templates/medoid_{number}_{label}/complex_heavy.pdb"
    for number, label in (
        ("01", "delivery_full_010_frame_068"),
        ("02", "delivery_full_011_frame_006"),
        ("03", "delivery_full_007_frame_095"),
    )
]


def write_new(path: Path, value: Any) -> None:
    """Write a new JSON artifact without replacing prior evidence."""
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def build_request(protection: dict[str, Any]) -> dict[str, Any]:
    """Return the frozen budgets and gates for this pilot, without execution."""
    mask = protection["mutable"]
    tied = symmetry_manifest(mask)
    return {
        "schema_version": 1,
        "status": "queued_pending_rfd3_and_execution_release",
        "campaign": "nyl12_relay_broad_pilot_20261126",
        "protection_receipt": "canonical_protection.json",
        "protection_sha256": None,
        "original_parent_ids": list(ORIGINAL_PARENTS),
        "original_parent_pdbs": {
            key: str(path) for key, path in ORIGINAL_PARENTS.items()
        },
        "remodeled_parent_gate": {
            "manifest": str(CAMPAIGN / "backbone_qualification/expanded_parent_manifest.json"),
            "selection_manifest": str(CAMPAIGN / "backbone_qualification/selected_parent_manifest.json"),
            "minimum_accepted": 0,
            "maximum_selected": 4,
            "require_exact_hashes_and_scientific_acceptance": True,
        },
        "generation": {
            "temperature": 0.1,
            "native_identity_bias": 0,
            "max_raw_samples": 24,
            "original_samples_per_parent": 6,
            "remodeled_samples_per_parent": 3,
            "batch_size": 1,
            "no_padding_or_retries": True,
            "allowed_amino_acids": "all standard amino acids at mutable sites",
            "mask": mask,
            **tied,
        },
        "execution_gate": {
            "accepted_parent_manifests_and_referenced_hashes_valid": False,
            "resource_release_present": False,
            "root_execution_authorization": True,
            "gpu": 0,
            "max_concurrent_jobs": 1,
        },
        "downstream": {
            "pool_complete_sequences_and_preserve_parent_lineage": True,
            "report_actual_mutation_counts_without_forced_caps": True,
            "rank_activity": False,
            "no_activity_or_energy_claim_from_likelihood": True,
        },
        "known_blockers": [
            "RFD3 qualification must pass for both initial parents before expansion.",
            "Current accepted_parent_manifest is not accepted; do not execute generation.",
            "Each RFD3 parent needs validated native LigandMPNN input, mapping, parser and PA6 graph receipts.",
        ],
    }


def prepare(output: Path) -> tuple[Path, Path]:
    """Write protection census and the non-executable bounded queue request."""
    output.mkdir(parents=True, exist_ok=False)
    protection = protection_census(CENSUS_PARENTS)
    protection["source_sha256"] = {
        str(path): sha256_file(path) for path in CENSUS_PARENTS
    }
    protection["method"] = {
        "rule": "fixed functional/maturation/disulfide positions plus recurring inter-repeat contacts",
        "interface_distance_angstrom": 4.0,
        "recurrence": "at least 3 of 4 equivalent copies in each of all 3 original TS parents",
        "excluded_intra_repeat_pairs": ["AE", "BF", "CG", "DH"],
        "no_rsa_filter": True,
        "interface_protection_is_operational_not_proven_essential": True,
    }
    protection_path = output / "canonical_protection.json"
    write_new(protection_path, protection)
    request = build_request(protection)
    request["protection_sha256"] = sha256_file(protection_path)
    request_path = output / "pilot_request.json"
    write_new(request_path, request)
    return protection_path, request_path


def main() -> None:
    """Write immutable preparation artifacts at a new output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    protection, request = prepare(args.output)
    print(f"protection={protection}")
    print(f"request={request}")


if __name__ == "__main__":
    main()
