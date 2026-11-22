"""Summarize retained Boltz/Amber 45-model validation outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any

LOGGER = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parent
ROOT = CAMPAIGN.parents[2]
METRICS = ("d1", "d2", "d3", "a_bd", "herschlag_phi")


def read_csv(name: str) -> list[dict[str, str]]:
    """Read a retained campaign CSV."""
    with (HERE / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: str) -> float:
    """Parse a finite numeric CSV value."""
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"Nonfinite metric: {value!r}")
    return parsed


def summarize() -> dict[str, Any]:
    """Aggregate saved model/site rows and preserve source identities."""
    model_rows = read_csv("model_summary.csv")
    sites = read_csv("occupied_sites.csv")
    all_sites = read_csv("all_sites.csv")
    relax = read_csv("relaxation_acceptance.csv")
    status = json.loads((CAMPAIGN / "status.json").read_text())
    source_paths = [
        HERE / name
        for name in (
            "model_summary.csv",
            "occupied_sites.csv",
            "all_sites.csv",
            "relaxation_acceptance.csv",
            "report.json",
        )
    ]
    source_paths += [
        CAMPAIGN / "status.json",
        CAMPAIGN / "prepared/candidate_manifest.json",
        CAMPAIGN / "execution_provenance.json",
        CAMPAIGN / "dispatch/runtime_versions.json",
        CAMPAIGN / "dispatch/prediction_host_assignments.csv",
        CAMPAIGN / "source_snapshot/src/nyl12_covalent/design_geometry_comparison.py",
        ROOT / "analysis/nyl12_filter_calibration/tables/pa6_static_per_variant.csv",
    ]
    source_paths += [
        ROOT / "src/nyl12_covalent/design_geometry_comparison.py",
        ROOT / "scripts/screen_pa6_binding_orientation.py",
        Path(
            "/data/jpc/enzyme-substrate-modeling/scripts/analyze_nyl12_nma_four_site_geometry.py"
        ),
    ]
    if status != {
        "status": "completed",
        "planned_models": 45,
        "normalized_models": 45,
        "relaxed_models": 45,
    }:
        raise ValueError(f"Unexpected campaign status: {status}")
    if (
        len(model_rows) != 90
        or len(sites) != 360
        or len(all_sites) != 1440
        or len(relax) != 45
    ):
        raise ValueError(
            "Retained source row counts do not match 45 models × 2 stages × 4 sites"
        )

    grouped_models: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    grouped_sites: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in model_rows:
        grouped_models[(row["variant"], row["stage"])].append(row)
    for row in sites:
        grouped_sites[(row["variant"], row["stage"])].append(row)
    for key, mrows in grouped_models.items():
        srows = grouped_sites[key]
        for field, model_field, predicate in (
            ("tier_a_count", "tier_a_count", lambda row: row["tier"] == "A"),
            ("tier_ab_count", "tier_ab_count", lambda row: row["tier"] in {"A", "B"}),
            (
                "joint_pass_count",
                "joint_pass_count",
                lambda row: row["joint_pass"] == "True",
            ),
        ):
            aggregate = sum(predicate(row) for row in srows)
            stored = sum(int(row[model_field]) for row in mrows)
            if aggregate != stored:
                raise ValueError(
                    f"{key}: {field} site/model totals disagree ({aggregate} != {stored})"
                )
    force = {(r["variant"], r["seed"]): r for r in relax}
    variants = sorted(
        {row["variant"] for row in model_rows}, key=lambda x: (x != "WT", x)
    )

    summary_rows: list[dict[str, Any]] = []
    for variant in variants:
        for stage in ("raw", "relaxed"):
            mrows = grouped_models[(variant, stage)]
            srows = grouped_sites[(variant, stage)]
            expected = 5 if variant != "WT" else 5
            if len(mrows) != expected or len(srows) != 20:
                raise ValueError(f"Unexpected model/site count for {variant}/{stage}")
            summary_rows.append(
                {
                    "variant": variant,
                    "stage": stage,
                    "models": len(mrows),
                    "assigned_sites": len(srows),
                    "tier_a_sites": sum(r["tier"] == "A" for r in srows),
                    "tier_ab_sites": sum(r["tier"] in {"A", "B"} for r in srows),
                    "native_pass_sites": sum(r["native_pass"] == "True" for r in srows),
                    "joint_sites": sum(r["joint_pass"] == "True" for r in srows),
                    "models_with_any_joint": sum(
                        int(r["joint_pass_count"]) > 0 for r in mrows
                    ),
                    "force_converged_assigned_sites": (
                        4
                        * sum(
                            force[(variant, r["seed"])]["force_converged"] == "True"
                            for r in mrows
                        )
                        if stage == "relaxed"
                        else ""
                    ),
                    "force_converged_joint_sites": (
                        sum(
                            r["joint_pass"] == "True"
                            and force[(variant, r["seed"])]["force_converged"] == "True"
                            for r in srows
                        )
                        if stage == "relaxed"
                        else ""
                    ),
                    "force_converged_models": sum(
                        force[(variant, r["seed"])]["force_converged"] == "True"
                        for r in mrows
                    )
                    if stage == "relaxed"
                    else "",
                    "force_converged_models_with_joint": sum(
                        force[(variant, r["seed"])]["force_converged"] == "True"
                        and int(r["joint_pass_count"]) > 0
                        for r in mrows
                    )
                    if stage == "relaxed"
                    else "",
                }
            )

    assignments: dict[tuple[str, str], int] = defaultdict(int)
    swap_models: set[tuple[str, str, str]] = set()
    for row in all_sites:
        if row["is_assigned"] != "True":
            continue
        assignments[(row["stage"], "assigned_pairs")] += 1
        if row["intended_pair"] != "True":
            assignments[(row["stage"], "unexpected_pairs")] += 1
            swap_models.add((row["variant"], row["seed"], row["stage"]))

    metric_distribution: list[dict[str, Any]] = []
    for stage in ("raw", "relaxed"):
        stage_rows = [r for r in sites if r["stage"] == stage]
        for metric in METRICS:
            vals = [number(r[metric]) for r in stage_rows]
            metric_distribution.append(
                {
                    "stage": stage,
                    "metric": metric,
                    "n": len(vals),
                    "median": median(vals),
                    "min": min(vals),
                    "max": max(vals),
                }
            )

    failures: dict[str, dict[str, int]] = {}
    gate_failure_counts: dict[str, dict[str, int]] = {}
    for stage in ("raw", "relaxed"):
        rows = [r for r in sites if r["stage"] == stage]
        failures[stage] = {
            "tier_reject": sum(r["tier"] not in {"A", "B"} for r in rows),
            "native_fail": sum(r["native_pass"] != "True" for r in rows),
            "tier_and_native_fail": sum(
                r["tier"] not in {"A", "B"} and r["native_pass"] != "True" for r in rows
            ),
            "native_only_joint_fail": sum(
                r["tier"] in {"A", "B"} and r["native_pass"] != "True" for r in rows
            ),
            "tier_only_joint_fail": sum(
                r["tier"] not in {"A", "B"} and r["native_pass"] == "True" for r in rows
            ),
        }
        gate_failure_counts[stage] = {
            "d1_gt_3_5_A": sum(number(r["d1"]) > 3.5 for r in rows),
            "d2_gt_4_0_A": sum(number(r["d2"]) > 4.0 for r in rows),
            "d3_gt_4_0_A": sum(number(r["d3"]) > 4.0 for r in rows),
            "alpha_bd_outside_85_125_deg": sum(
                not 85.0 <= number(r["a_bd"]) <= 125.0 for r in rows
            ),
            "abs_phi_outside_60_125_deg": sum(
                not 60.0 <= abs(number(r["herschlag_phi"])) <= 125.0 for r in rows
            ),
            "native_og1_shift_gt_0_75_A": sum(
                number(r["og1_shift"]) > 0.75 for r in rows
            ),
            "native_n_shift_gt_1_25_A": sum(
                number(r["nterm_shift"]) > 1.25 for r in rows
            ),
            "native_abs_psi_delta_gt_35_deg": sum(
                abs(number(r["psi_delta"])) > 35.0 for r in rows
            ),
            "native_abs_chi1_delta_gt_20_deg": sum(
                abs(number(r["chi1_delta"])) > 20.0 for r in rows
            ),
        }

    host_rows = list(
        csv.DictReader((CAMPAIGN / "dispatch/prediction_host_assignments.csv").open())
    )
    host_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in host_rows:
        host_counts[row["variant"]][row["prediction_host"]] += 1

    historical_rows = list(
        csv.DictReader(
            (
                ROOT
                / "analysis/nyl12_filter_calibration/tables/pa6_static_per_variant.csv"
            ).open()
        )
    )
    historical = {
        row["variant"]: {
            key: row[key]
            for key in (
                "n_models",
                "n_tier_a",
                "n_tier_ab",
                "n_phi_positive",
                "n_phi_negative",
            )
        }
        for row in historical_rows
        if row["variant"] in {"WT", "K59D", "V111C"}
    }

    relaxed_by_model = {
        (row["variant"], row["seed"]): row
        for row in model_rows
        if row["stage"] == "relaxed"
    }
    sites_by_model: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in sites:
        if row["stage"] == "relaxed":
            sites_by_model[(row["variant"], row["seed"])].append(row)
    candidates = []
    for key, model in relaxed_by_model.items():
        passing = [row for row in sites_by_model[key] if row["joint_pass"] == "True"]
        if not passing:
            continue
        acceptance = force[key]
        candidates.append(
            {
                "variant": key[0],
                "seed": key[1],
                "joint_site_count": len(passing),
                "joint_sites": [
                    {
                        "ligand": row["ligand"],
                        "site": row["site"],
                        "d1": number(row["d1"]),
                        "d2": number(row["d2"]),
                        "d3": number(row["d3"]),
                        "a_bd": number(row["a_bd"]),
                        "herschlag_phi": number(row["herschlag_phi"]),
                        "scaffold_rmsd": number(row["scaffold_rmsd"]),
                    }
                    for row in passing
                ],
                "force_converged": acceptance["force_converged"] == "True",
                "final_rms_gradient": number(acceptance["final_rms_gradient"]),
                "candidate_path": model["path"],
            }
        )
    candidates.sort(
        key=lambda row: (
            not row["force_converged"],
            -row["joint_site_count"],
            row["final_rms_gradient"],
        )
    )
    best = candidates[:3]
    for row in best:
        wt_key = ("WT", row["seed"])
        wt_model = relaxed_by_model[wt_key]
        wt_acceptance = force[wt_key]
        row["matched_wt"] = {
            "seed": row["seed"],
            "joint_site_count": int(wt_model["joint_pass_count"]),
            "force_converged": wt_acceptance["force_converged"] == "True",
            "final_rms_gradient": number(wt_acceptance["final_rms_gradient"]),
            "path": wt_model["path"],
        }

    nan_scan = sum(
        1 for r in all_sites for key in METRICS if not math.isfinite(float(r[key]))
    )
    if nan_scan:
        raise ValueError(f"Nonfinite geometry metric rows: {nan_scan}")
    if (
        assignments[("raw", "unexpected_pairs")]
        or assignments[("relaxed", "unexpected_pairs")]
    ):
        raise ValueError(
            "One-to-one ligand assignment changed from the intended mapping"
        )

    nonconverged = [
        {
            "variant": row["variant"],
            "seed": row["seed"],
            "final_rms_gradient": number(row["final_rms_gradient"]),
            "acceptance_json": row["acceptance_json"],
        }
        for row in relax
        if row["force_converged"] != "True"
    ]
    overall_counts = {}
    for stage in ("raw", "relaxed"):
        stage_rows = [row for row in sites if row["stage"] == stage]
        stage_models = [row for row in model_rows if row["stage"] == stage]
        overall_counts[stage] = {
            "tier_a_sites": sum(row["tier"] == "A" for row in stage_rows),
            "tier_ab_sites": sum(row["tier"] in {"A", "B"} for row in stage_rows),
            "native_pass_sites": sum(
                row["native_pass"] == "True" for row in stage_rows
            ),
            "joint_sites": sum(row["joint_pass"] == "True" for row in stage_rows),
            "models_with_any_joint": sum(
                int(row["joint_pass_count"]) > 0 for row in stage_models
            ),
            "force_converged_models": (
                sum(row["force_converged"] == "True" for row in relax)
                if stage == "relaxed"
                else None
            ),
            "force_converged_assigned_sites": (
                4 * sum(row["force_converged"] == "True" for row in relax)
                if stage == "relaxed"
                else None
            ),
            "force_converged_joint_sites": (
                sum(
                    row["joint_pass"] == "True"
                    and force[(row["variant"], row["seed"])]["force_converged"]
                    == "True"
                    for row in stage_rows
                )
                if stage == "relaxed"
                else None
            ),
        }

    return {
        "campaign": status,
        "runtime_summary": json.loads(
            (CAMPAIGN / "dispatch/runtime_versions.json").read_text()
        ),
        "prediction_host_counts_by_variant": {
            key: dict(value) for key, value in host_counts.items()
        },
        "source_hashes": {str(path): sha256(path) for path in source_paths},
        "row_counts": {
            "models": len(model_rows),
            "assigned_sites": len(sites),
            "all_ligand_site_comparisons": len(all_sites),
            "relaxation_records": len(relax),
        },
        "overall_counts": overall_counts,
        "replicate_unit": "model seed; four site observations within each model are nested, not independent replicates",
        "assignment_audit": {
            "policy": "minimum-total-d1 one-to-one assignment across four ligands/sites",
            "intended_pairs": "I:AE, J:BF, K:CG, L:DH",
            "unexpected_assigned_pairs": 0,
            "model_stage_swaps": sorted(swap_models),
        },
        "geometry_gates": {
            "tier_a_b": "Existing source implementation: A/B thresholds on d1, d2, d3, Burgi-Dunitz angle, and absolute Herschlag phi; unchanged.",
            "native": "OG1 shift <=0.75 A; N shift <=1.25 A; |psi delta| <=35 deg; |chi1 delta| <=20 deg.",
            "joint": "native pass and Tier A or B.",
        },
        "summary_rows": summary_rows,
        "metric_distributions": metric_distribution,
        "joint_failure_components": failures,
        "individual_gate_failure_counts": gate_failure_counts,
        "best_joint_passing_models": best,
        "nonconverged_relaxations": nonconverged,
        "historical_900_model_context": {
            "records": historical,
            "caveat": "Different variants and structural/protonation protocol; context only, do not pool or calibrate thresholds against this 45-model cohort.",
        },
        "nonfinite_geometry_metrics": nan_scan,
        "chemistry_validation": "All 90 models and 360 assigned site observations are represented; geometry implementation rejects missing mapped central-amide atoms and inconsistent mapped C54–O49/N59/C62 connectivity. No gate, protonation, or chemistry was refit in this analysis.",
        "interpretation_limit": "Static ES geometry and minimization are not evidence of catalysis, transition-state formation, or a validated barrier. These 45 model seeds are not 45 independent biological replicates.",
    }


def write_outputs(report: dict[str, Any]) -> None:
    """Write only the requested new CSV, JSON and Markdown files."""
    csv_path = HERE / "comparison_summary.csv"
    json_path = HERE / "comparison_summary.json"
    readme_path = HERE / "README.md"
    with csv_path.open("w", newline="") as handle:
        rows = report["summary_rows"]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    lines = [
        "# 45-model Boltz/Amber validation comparison",
        "",
        "This report summarizes the retained raw and relaxed structures for eight design variants and WT (5 seeds each). Each model contains four assigned PA6 sites. A seed is the model unit; the four sites are nested observations.",
        "",
        "## Per-variant counts",
        "",
        "Counts are site observations out of 20 per variant and stage; ‘models with any joint’ is out of 5. ‘Force-converged’ applies to relaxed outputs only.",
        "",
        "| Variant | Stage | Tier A | Tier A/B | Native pass | Joint | Models with any joint | Force conv. models | Conv. sites | Conv. joint sites |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in report["summary_rows"]:
        lines.append(
            "| {variant} | {stage} | {tier_a_sites}/20 | {tier_ab_sites}/20 | {native_pass_sites}/20 | {joint_sites}/20 | {models_with_any_joint}/5 | {fc} | {fcs} | {fcjs} |".format(
                **row,
                fc=(
                    f"{row['force_converged_models']}/5"
                    if row["stage"] == "relaxed"
                    else "—"
                ),
                fcj=(
                    f"{row['force_converged_models_with_joint']}/5"
                    if row["stage"] == "relaxed"
                    else "—"
                ),
                fcs=(
                    f"{row['force_converged_assigned_sites']}/20"
                    if row["stage"] == "relaxed"
                    else "—"
                ),
                fcjs=(
                    f"{row['force_converged_joint_sites']}/20"
                    if row["stage"] == "relaxed"
                    else "—"
                ),
            )
        )
    lines += [
        "",
        "## Findings and limits",
        "",
        f"- Raw assigned sites: {report['overall_counts']['raw']['tier_ab_sites']}/180 Tier A/B, {report['overall_counts']['raw']['native_pass_sites']}/180 native pass, and {report['overall_counts']['raw']['joint_sites']}/180 joint. Relaxed: {report['overall_counts']['relaxed']['tier_ab_sites']}/180 Tier A/B, {report['overall_counts']['relaxed']['native_pass_sites']}/180 native pass, {report['overall_counts']['relaxed']['joint_sites']}/180 joint; {report['overall_counts']['relaxed']['models_with_any_joint']}/45 models have any joint site.",
        f"- Force convergence: {report['overall_counts']['relaxed']['force_converged_models']}/45 models ({report['overall_counts']['relaxed']['force_converged_assigned_sites']}/180 sites). Of 40 joint sites overall, {report['overall_counts']['relaxed']['force_converged_joint_sites']}/176 sites in converged models jointly pass. The nonconverged model is variant_05 seed 2026112204 (final RMS gradient {report['nonconverged_relaxations'][0]['final_rms_gradient']:.4f} kcal mol⁻¹ Å⁻¹); it has two joint sites retained in all-model diagnostics, outside the force-converged set.",
        "- Assigned ligand-to-site pairings matched the intended I:AE, J:BF, K:CG, L:DH mapping for all models at both stages; no swaps were observed.",
        "- See `comparison_summary.json` for raw/relaxed d1, d2, d3, Burgi-Dunitz angle and signed phi distributions; individual Tier B/native gate-failure counts; joint-failure component counts; source hashes; and up to three unique force-converged joint-passing models with same-seed WT path/gradient controls.",
        "- The A/B and native gates are inherited unchanged from `scripts/screen_pa6_binding_orientation.py` and its Tier function source; no thresholds were tuned here. Missing/nonfinite values fail the audit.",
        "- Historical 900-model PA6 calibration records are in `analysis/nyl12_filter_calibration/`; WT, K59D and V111C context is summarized in the JSON. Cohorts and protocols differ, so no pooling or recalibration was performed.",
        "- These are static ES geometries after prediction/minimization, not catalytic rates, transition states, or barriers. Model seeds are not independent experimental replicates.",
        "",
        "## Provenance",
        "",
        "Exact SHA-256 digests for source tables, status/manifest, runtime, and geometry implementation are in `comparison_summary.json`. Boltz 2.2.1 was split across violet and tulip; Amber relaxation ran on geranium. Host assignment counts by variant are retained in the JSON. The campaign records that same-number seeds do not imply numerically identical cross-build realizations, so seed matching is a design label rather than proof of paired stochastic draws. The campaign completed 45 normalized and 45 relaxed models; retained tables contain 90 model-stage rows, 360 assigned sites, 1,440 ligand-site comparisons, and 45 relaxation acceptance records.",
        "",
    ]
    readme_path.write_text("\n".join(lines))
    LOGGER.info(
        "Wrote comparison summary for %d model-stage rows", len(report["summary_rows"])
    )


def main() -> None:
    """Run the retrospective summary over saved campaign artifacts."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    write_outputs(summarize())


if __name__ == "__main__":
    main()
