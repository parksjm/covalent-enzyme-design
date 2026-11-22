"""Summarize and retain provenance-bound structural-diversity results."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from diversity_geometry import load_models, pairwise_rmsds, sha256, verify_fit_helper

LOGGER = logging.getLogger(__name__)
PAIR_FIELDS = """variant stage seed_left seed_right force_converged_left
force_converged_right both_force_converged global_ca_rmsd_angstrom
active_ca_rmsd_angstrom""".split()
SUMMARY_FIELDS = """variant stage pair_count global_median global_maximum
active_median active_maximum force_converged_pair_count
force_converged_global_median force_converged_global_maximum
force_converged_active_median force_converged_active_maximum
delta_global_median_vs_WT delta_global_maximum_vs_WT
delta_active_median_vs_WT delta_active_maximum_vs_WT joint_sites_of_20
force_converged_joint_sites models_with_any_joint_of_5""".split()


def _distribution(values: list[float]) -> dict[str, float | int]:
    if not values or not np.isfinite(values).all():
        raise ValueError("RMSD distribution is empty or nonfinite")
    return {
        "pair_count": len(values),
        "median": float(np.median(values)),
        "maximum": float(max(values)),
    }


def _joint_pass_context(path: Path) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["stage"] == "relaxed":
                result[row["variant"]] = {
                    "joint_sites_of_20": int(row["joint_sites"]),
                    "force_converged_joint_sites": int(
                        row["force_converged_joint_sites"]
                    ),
                    "models_with_any_joint_of_5": int(row["models_with_any_joint"]),
                }
    return result


def _summaries(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for variant in sorted({row["variant"] for row in rows}):
        for stage in ("raw", "relaxed"):
            group = [
                row
                for row in rows
                if row["variant"] == variant and row["stage"] == stage
            ]
            summary: dict[str, Any] = {
                "variant": variant,
                "stage": stage,
                "global_ca_rmsd_angstrom": _distribution(
                    [row["global_ca_rmsd_angstrom"] for row in group]
                ),
                "active_ca_rmsd_angstrom": _distribution(
                    [row["active_ca_rmsd_angstrom"] for row in group]
                ),
            }
            if stage == "relaxed":
                converged = [row for row in group if row["both_force_converged"]]
                summary["both_force_converged_pairs"] = {
                    "global_ca_rmsd_angstrom": _distribution(
                        [row["global_ca_rmsd_angstrom"] for row in converged]
                    ),
                    "active_ca_rmsd_angstrom": _distribution(
                        [row["active_ca_rmsd_angstrom"] for row in converged]
                    ),
                }
            summaries.append(summary)
    metrics = ("global_ca_rmsd_angstrom", "active_ca_rmsd_angstrom")
    for stage in ("raw", "relaxed"):
        wt = next(
            row for row in summaries if row["variant"] == "WT" and row["stage"] == stage
        )
        for row in summaries:
            if row["stage"] != stage or row["variant"] == "WT":
                continue
            row["delta_vs_WT"] = {
                metric: {
                    stat: row[metric][stat] - wt[metric][stat]
                    for stat in ("median", "maximum")
                }
                for metric in metrics
            }
            if stage == "relaxed":
                row["both_force_converged_pairs"]["delta_vs_WT"] = {
                    metric: {
                        stat: row["both_force_converged_pairs"][metric][stat]
                        - wt["both_force_converged_pairs"][metric][stat]
                        for stat in ("median", "maximum")
                    }
                    for metric in metrics
                }
    return summaries


def run_analysis(
    root: Path, boltz_root: Path, output_dir: Path, script_path: Path
) -> dict[str, Any]:
    """Analyze 45 retained models and write compact JSON/CSV summaries.

    Args:
        root: Repository root.
        boltz_root: Completed 45-model validation root.
        output_dir: New or existing analysis directory containing this script.
        script_path: CLI script whose hash is retained as executable provenance.

    Returns:
        JSON-ready analysis report with per-model, pairwise, and summary records.

    Raises:
        FileExistsError: If a result file already exists.
        FileNotFoundError: If any expected structure or receipt is absent.
        ValueError: If model lineage, sequences, alignment, or measurements fail.
    """
    outputs = [
        output_dir / "diversity_summary.json",
        output_dir / "pairwise_rmsd.csv",
        output_dir / "variant_summary.csv",
    ]
    existing = next((path for path in outputs if path.exists()), None)
    if existing:
        raise FileExistsError(existing)
    helper = root / "src/nyl12_covalent/catalytic_alignment_geometry.py"
    geometry_source = output_dir / "diversity_geometry.py"
    reporting_source = output_dir / "diversity_reporting.py"
    mask_path = (
        root
        / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/preflight/parent_and_mask_proposal.json"
    )
    summary_path = boltz_root / "analysis/comparison_summary.csv"
    models, hashes = load_models(boltz_root)
    sources = [
        boltz_root / "prepared/candidate_manifest.json",
        boltz_root / "status.json",
        boltz_root / "completion.json",
        summary_path,
        mask_path,
        helper,
        geometry_source,
        reporting_source,
        script_path.resolve(),
    ]
    hashes.update({str(path.resolve()): sha256(path) for path in sources})
    pairs = pairwise_rmsds(models)
    summaries = _summaries(pairs)
    context = _joint_pass_context(summary_path)
    for summary in summaries:
        if summary["stage"] == "relaxed":
            summary["structural_joint_pass_context"] = context[summary["variant"]]
    model_status = [
        {
            "variant": model["variant"],
            "seed": model["seed"],
            "force_converged": model["force_converged"],
            "final_rms_gradient": model["final_rms_gradient"],
            "raw_path": str(model["raw_path"].resolve()),
            "raw_sha256": model["raw_sha256"],
            "relaxed_path": str(model["relaxed_path"].resolve()),
            "relaxed_sha256": model["relaxed_sha256"],
        }
        for model in models
    ]
    report: dict[str, Any] = {
        "analysis_definition": {
            "global_fit": "proper-rotation Kabsch fit on all 1,312 A-H protein CAs in canonical chain/position order",
            "active_subset": {
                "alpha_positions_per_A_to_D_copy": [
                    50,
                    51,
                    52,
                    73,
                    76,
                    101,
                    111,
                    112,
                    113,
                    114,
                    115,
                    116,
                    151,
                    152,
                    153,
                    156,
                    185,
                    186,
                    187,
                    188,
                ],
                "beta_positions_per_E_to_H_copy": [1, 39, 40, 41, 42],
                "ca_atom_count": 100,
                "active_rmsd_uses_global_fit": True,
            },
            "pairs": "all ten unordered within-variant seed pairs for each stage; pairwise values are not independent model replicates",
            "relaxed_outputs": "all completed models included; force-converged-only statistics use pairs where both models meet the recorded force tolerance",
            "scope": "descriptive static structural variation; no catalytic activity or causal heterogeneity inference",
        },
        "fit_validation": verify_fit_helper(),
        "model_count": len(model_status),
        "pairwise_record_count": len(pairs),
        "force_converged_models": sum(
            model["force_converged"] for model in model_status
        ),
        "nonconverged_models": [
            model for model in model_status if not model["force_converged"]
        ],
        "models": model_status,
        "variant_summaries": summaries,
        "source_sha256": dict(sorted(hashes.items())),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs[0].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    with outputs[1].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_FIELDS)
        writer.writeheader()
        writer.writerows(pairs)
    with outputs[2].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for summary in summaries:
            global_stats = summary["global_ca_rmsd_angstrom"]
            active_stats = summary["active_ca_rmsd_angstrom"]
            converged = summary.get("both_force_converged_pairs", {})
            delta = summary.get("delta_vs_WT", {})
            writer.writerow(
                {
                    "variant": summary["variant"],
                    "stage": summary["stage"],
                    "pair_count": global_stats["pair_count"],
                    "global_median": global_stats["median"],
                    "global_maximum": global_stats["maximum"],
                    "active_median": active_stats["median"],
                    "active_maximum": active_stats["maximum"],
                    "force_converged_pair_count": converged.get(
                        "global_ca_rmsd_angstrom", {}
                    ).get("pair_count", ""),
                    "force_converged_global_median": converged.get(
                        "global_ca_rmsd_angstrom", {}
                    ).get("median", ""),
                    "force_converged_global_maximum": converged.get(
                        "global_ca_rmsd_angstrom", {}
                    ).get("maximum", ""),
                    "force_converged_active_median": converged.get(
                        "active_ca_rmsd_angstrom", {}
                    ).get("median", ""),
                    "force_converged_active_maximum": converged.get(
                        "active_ca_rmsd_angstrom", {}
                    ).get("maximum", ""),
                    "delta_global_median_vs_WT": delta.get(
                        "global_ca_rmsd_angstrom", {}
                    ).get("median", ""),
                    "delta_global_maximum_vs_WT": delta.get(
                        "global_ca_rmsd_angstrom", {}
                    ).get("maximum", ""),
                    "delta_active_median_vs_WT": delta.get(
                        "active_ca_rmsd_angstrom", {}
                    ).get("median", ""),
                    "delta_active_maximum_vs_WT": delta.get(
                        "active_ca_rmsd_angstrom", {}
                    ).get("maximum", ""),
                    **summary.get("structural_joint_pass_context", {}),
                }
            )
    LOGGER.info(
        "Wrote %d pairwise records for %d models", len(pairs), len(model_status)
    )
    return report
