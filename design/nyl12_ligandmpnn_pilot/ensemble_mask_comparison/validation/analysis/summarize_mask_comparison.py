#!/usr/bin/env python3
"""Build an auditable summary of the completed mask comparison."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[5]
VALIDATION = (
    ROOT
    / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation"
)
ANALYSIS = VALIDATION / "analysis"
MANIFEST_PATH = VALIDATION / "prepared/candidate_manifest.json"
EXPECTED_ASSIGNMENT = {"I:1": "AE", "J:1": "BF", "K:1": "CG", "L:1": "DH"}
EXPECTED_SITES = {"AE", "BF", "CG", "DH"}
SEEDS = {2026112201, 2026112202, 2026112203, 2026112204, 2026112205}


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV file into dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def truth(value: str) -> bool:
    """Parse a serialized boolean without accepting ambiguous values."""
    if value == "True":
        return True
    if value == "False":
        return False
    raise ValueError(f"Expected True or False, got {value!r}")


def candidate_groups(manifest: dict[str, Any]) -> dict[str, str]:
    """Map each candidate to its comparison arm or control cohort."""
    result: dict[str, str] = {}
    for candidate in manifest["candidates"]:
        arms = candidate["arms"]
        if len(arms) > 1:
            raise ValueError(f"Candidate has multiple arms: {candidate['id']}")
        result[candidate["id"]] = arms[0] if arms else candidate["cohort_tag"]
    return result


def summarize_sites(rows: list[dict[str, str]]) -> dict[str, int]:
    """Count geometry gates and signed-phi outcomes for intended sites."""
    joint = [row for row in rows if row["tier"] in {"A", "B"} and truth(row["native_pass"])]
    if any(truth(row["joint_pass"]) != (row in joint) for row in rows):
        raise ValueError("Saved joint_pass disagrees with Tier A/B and native_pass")
    by_model: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_model[(row["variant"], row["seed"])].append(row)
    joint_ids = {id(row) for row in joint}
    signs = {
        key: {row["orientation"] for row in model_rows if id(row) in joint_ids}
        for key, model_rows in by_model.items()
    }
    return {
        "candidate_count": len({row["variant"] for row in rows}),
        "model_count": len(by_model),
        "site_count": len(rows),
        "tier_a_sites": sum(row["tier"] == "A" for row in rows),
        "tier_ab_sites": sum(row["tier"] in {"A", "B"} for row in rows),
        "native_pass_sites": sum(truth(row["native_pass"]) for row in rows),
        "joint_pass_sites": len(joint),
        "models_any_joint_pass": sum(bool(value) for value in signs.values()),
        "negative_phi_joint_sites": sum(
            row["orientation"] == "negative_phi" for row in joint
        ),
        "positive_phi_joint_sites": sum(
            row["orientation"] == "positive_phi" for row in joint
        ),
        "negative_only_joint_models": sum(value == {"negative_phi"} for value in signs.values()),
    }


def validate_sources(
    manifest: dict[str, Any],
    occupied: list[dict[str, str]],
    model_rows: list[dict[str, str]],
) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    """Validate expected models, assignments, source counts, and relaxations."""
    groups = candidate_groups(manifest)
    designs = [c for c in manifest["candidates"] if c["arms"]]
    controls = [c for c in manifest["candidates"] if not c["arms"] and c["id"] != "WT"]
    if len(designs) != 8 or len(controls) != 2 or groups.get("WT") != "baseline":
        raise ValueError("Expected eight designs, two controls, and WT baseline")
    if sum(groups[c["id"]] == "broader_29_positions" for c in designs) != 4:
        raise ValueError("Expected four broader-arm candidates")
    if sum(groups[c["id"]] == "focused_19_positions" for c in designs) != 4:
        raise ValueError("Expected four focused-arm candidates")

    relaxed_models = [row for row in model_rows if row["stage"] == "relaxed"]
    if len(relaxed_models) != 55:
        raise ValueError(f"Expected 55 relaxed models, found {len(relaxed_models)}")
    expected_ids = {candidate["id"] for candidate in manifest["candidates"]}
    if {row["variant"] for row in relaxed_models} != expected_ids:
        raise ValueError("Relaxed model candidate IDs differ from manifest")
    by_model: dict[tuple[str, str], dict[str, str]] = {}
    for model in relaxed_models:
        key = (model["variant"], model["seed"])
        if key in by_model:
            raise ValueError(f"Duplicate relaxed model summary: {key}")
        assignment = json.loads(model["assignment"])
        if assignment != EXPECTED_ASSIGNMENT:
            raise ValueError(f"Unexpected assignment for {key}: {assignment}")
        by_model[key] = model
    sites_by_model: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in occupied:
        if row["stage"] != "relaxed" or not truth(row["intended_pair"]):
            continue
        key = (row["variant"], row["seed"])
        if row["site"] != EXPECTED_ASSIGNMENT[row["ligand"]]:
            raise ValueError(f"Unexpected intended site for {key}: {row['site']}")
        if row["assigned_site"] != row["site"] or not truth(row["is_assigned"]):
            raise ValueError(f"Intended site was not assigned for {key}: {row}")
        sites_by_model[key].append(row)

    reused = {
        (item["variant"], str(item["seed"])): Path(item["acceptance_json"])
        for item in manifest["reused_models"]
    }
    acceptance_hashes: dict[str, str] = {}
    for key, model in by_model.items():
        seed = int(key[1])
        if seed not in SEEDS:
            raise ValueError(f"Unexpected seed for {key}")
        site_rows = sites_by_model[key]
        if len(site_rows) != 4 or {row["site"] for row in site_rows} != EXPECTED_SITES:
            raise ValueError(f"Expected four unique assigned sites for {key}")
        counts = summarize_sites(site_rows)
        if counts["site_count"] != int(model["ligand_count"]):
            raise ValueError(f"Ligand/site count mismatch for {key}")
        if counts["tier_a_sites"] != int(model["tier_a_count"]):
            raise ValueError(f"Tier A count mismatch for {key}")
        if counts["tier_ab_sites"] != int(model["tier_ab_count"]):
            raise ValueError(f"Tier A/B count mismatch for {key}")
        if counts["joint_pass_sites"] != int(model["joint_pass_count"]):
            raise ValueError(f"Joint count mismatch for {key}")
        acceptance = reused.get(key, Path(model["path"]).parent / "acceptance.json")
        if not acceptance.is_file():
            raise FileNotFoundError(f"Missing acceptance record for {key}: {acceptance}")
        record = json.loads(acceptance.read_text(encoding="utf-8"))
        if record.get("force_converged") is not True or record.get("process_exit_code") != 0:
            raise ValueError(f"Relaxation acceptance failed for {key}: {acceptance}")
        acceptance_hashes[str(acceptance)] = sha256(acceptance)

    if any(
        {int(seed) for variant, seed in by_model if variant == candidate}
        != SEEDS
        for candidate in expected_ids
    ):
        raise ValueError("At least one candidate is missing a required seed")
    if manifest["new_prediction_job_count"] != 52 or len(reused) != 3:
        raise ValueError("Expected 52 new predictions and three reused WT models")
    if len(acceptance_hashes) != 55:
        raise ValueError(f"Expected 55 distinct acceptance records, found {len(acceptance_hashes)}")
    return groups, acceptance_hashes


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write summary rows as a deterministic CSV table."""
    if not rows:
        raise ValueError(f"Refusing to write empty summary: {path}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Validate source artifacts and write compact summary reports."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    occupied_path = ANALYSIS / "occupied_sites.csv"
    model_path = ANALYSIS / "model_summary.csv"
    occupied = read_csv(occupied_path)
    model_rows = read_csv(model_path)
    groups, acceptance_hashes = validate_sources(manifest, occupied, model_rows)
    intended = [row for row in occupied if row["stage"] == "relaxed" and truth(row["intended_pair"])]
    candidate_rows: list[dict[str, Any]] = []
    for variant in sorted(groups):
        rows = [row for row in intended if row["variant"] == variant]
        summary = summarize_sites(rows)
        summary.update({"candidate": variant, "group": groups[variant]})
        candidate_rows.append(summary)
    if next(row for row in candidate_rows if row["candidate"] == "control_K59D")[
        "joint_pass_sites"
    ] != 0:
        raise ValueError("Known active control_K59D no longer has zero joint-pass sites")
    group_keys = ["broader_29_positions", "focused_19_positions", "baseline", "matched_positive_control"]
    group_rows: list[dict[str, Any]] = []
    for group in group_keys:
        rows = [row for row in intended if groups[row["variant"]] == group]
        summary = summarize_sites(rows)
        summary.update({"group": group})
        group_rows.append(summary)
    source_paths = [MANIFEST_PATH, occupied_path, model_path]
    source_hashes = {str(path): sha256(path) for path in source_paths}
    report = {
        "schema_version": 1,
        "campaign": "nyl12_ligandmpnn_pilot/ensemble_mask_comparison",
        "scope": "52 new predictions plus three reused WT models; distinct from the earlier 45-model campaign",
        "replicate_unit": "Prediction seed; sites within a model are nested, not independent replicates",
        "interpretation": {
            "joint_pass": "Tier A or B and native_pass",
            "small_ensemble": "Five seeds per candidate; this comparison does not support activity claims",
            "gate_sensitivity": "Known active control_K59D has zero joint-pass sites in this screen",
            "phi": "Opposite signed phi values indicate different local pose orientations; retain signs",
        },
        "counts": {
            "new_predictions": 52,
            "reused_wt_models": 3,
            "total_relaxed_models": 55,
            "design_candidates": 8,
            "matched_controls": 2,
            "wt_models": 5,
            "seeds": sorted(SEEDS),
            "sites_per_model": 4,
            "intended_assignment": EXPECTED_ASSIGNMENT,
            "all_force_converged_and_process_exit_zero": True,
        },
        "groups": group_rows,
        "candidates": candidate_rows,
        "source_sha256": source_hashes,
        "acceptance_sha256": acceptance_hashes,
    }
    (ANALYSIS / "comparison_summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(ANALYSIS / "comparison_summary.csv", candidate_rows)
    write_csv(ANALYSIS / "group_summary.csv", group_rows)
    logging.info("Wrote candidate and group summaries after validating 55 relaxations")


if __name__ == "__main__":
    main()
