"""Publish matched-cohort all-site assignments and complete amide geometry."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
from typing import Any

from nyl12_covalent.negative_amide_gate import evaluate_signed_amide_gate

SITES = ("AE", "BF", "CG", "DH")
METRICS = (
    "ligand",
    "carbonyl_c",
    "carbonyl_o",
    "leaving_n",
    "tier",
    "d1",
    "d2",
    "d3",
    "a_bd",
    "herschlag_phi",
)


def _read(path: Path, pins: dict[str, str]) -> dict[str, Any]:
    content = path.read_bytes()
    pins[str(path.resolve())] = hashlib.sha256(content).hexdigest()
    return json.loads(content)


def _table(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    if not columns:
        columns = ["case", "snapshot", "site", "measurement_status"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value) if isinstance(value, (list, dict)) else value
                    for key, value in row.items()
                }
            )


def _site_rows(
    case: str, snapshot: str, report: dict[str, Any]
) -> list[dict[str, Any]]:
    rows = report["all_amide_site_rows"]
    if report["target_phi_sign"] != 1:
        raise ValueError(f"Positive target required: {case}/{snapshot}")
    gate = evaluate_signed_amide_gate(rows, target_phi_sign=1)
    # Existing one-to-one assignment is authoritative; independent recomputation
    # detects incomplete or altered reports, without introducing another filter.
    if any(
        report[key] != gate[key]
        for key in ("assigned_site_order", "assigned", "assigned_ab_count")
    ):
        raise ValueError(f"Saved gate assignment inconsistent: {case}/{snapshot}")
    assigned = dict(zip(gate["assigned_site_order"], gate["assigned"], strict=True))
    result = []
    for site in SITES:
        selected = assigned[site]
        diagnostic = min(
            (r for r in rows if r["site"] == site),
            key=lambda r: (float(r["d1"]), str(r["ligand"]), str(r["carbonyl_c"])),
        )
        eligible_count = sum(r["site"] == site for r in gate["eligible_ab"])
        row = {
            "case": case,
            "snapshot": snapshot,
            "site": site,
            "measurement_status": "measured",
            "positive_ab_assigned": selected is not None,
            "assignment_status": "passing_positive_ab"
            if selected
            else "reject_no_passing_amide"
            if not eligible_count
            else "unassigned_one_to_one_conflict",
            "eligible_positive_ab_count": eligible_count,
            "selection_policy": gate["policy"],
        }
        row.update(
            {
                f"selected_{key}": selected.get(key) if selected else None
                for key in METRICS
            }
        )
        row.update(
            {f"diagnostic_nearest_d1_{key}": diagnostic.get(key) for key in METRICS}
        )
        result.append(row)
    return result


def summarize_cohort(cohort_root: Path, output: Path) -> dict[str, Any]:
    """Write new cohort tables without altering native or Amber evidence.

    Args:
        cohort_root: Explicit campaign root containing controller/summary.json.
        output: New comparison directory; existing outputs are never overwritten.

    Returns:
        Case outcomes, measurement counts, signed assignments and input hashes.

    Raises:
        ValueError: Complete cases lack valid all48 positive-phi gate evidence.
        FileExistsError: The requested output directory already exists.
    """
    if output.exists():
        raise FileExistsError(output)
    pins: dict[str, str] = {}
    controller = _read(cohort_root / "controller/summary.json", pins)
    site_rows, all_rows, cases = [], [], {}
    for case, record in sorted(controller["cases"].items()):
        if not case or Path(case).name != case or case in {".", ".."}:
            raise ValueError("Invalid case identity")
        case_summary: dict[str, Any] = {
            "status": record["status"],
            "mm_result_available": record["status"] == "completed",
        }
        if record["status"] != "completed":
            # Missing results are explicit placeholders, never four valid rejects.
            case_summary["measurement_status"] = "no_mm_result"
            for snapshot in ("before", "after"):
                site_rows.extend(
                    {
                        "case": case,
                        "snapshot": snapshot,
                        "site": site,
                        "measurement_status": "no_mm_result",
                        "assignment_status": "not_evaluated",
                        "positive_ab_assigned": None,
                        "case_status": record["status"],
                    }
                    for site in SITES
                )
        else:
            postprocess = Path(record["postprocess"])
            for snapshot in ("before", "after"):
                report = _read(postprocess / f"pose_gate_{snapshot}.json", pins)
                if report["case"] != case or report["snapshot"] != snapshot:
                    raise ValueError("Gate case/snapshot identity mismatch")
                measured = _site_rows(case, snapshot, report)
                site_rows.extend(measured)
                case_summary[f"positive_ab_{snapshot}"] = sum(
                    r["positive_ab_assigned"] for r in measured
                )
                all_rows.extend(
                    {**row, "case": case, "snapshot": snapshot}
                    for row in report["all_amide_site_rows"]
                )
        cases[case] = case_summary
    summary = {
        "target_phi_sign": 1,
        "cases": cases,
        "measured_mm_cases": sum(r["mm_result_available"] for r in cases.values()),
        "cases_without_mm_result": sum(
            not r["mm_result_available"] for r in cases.values()
        ),
        "measured_site_snapshot_rows": sum(
            r["measurement_status"] == "measured" for r in site_rows
        ),
        "all_amide_geometry_rows": len(all_rows),
        "input_sha256": pins,
        "native_thr_policy": "diagnostic only; no acceptance threshold added",
    }
    output.mkdir(parents=True)
    _table(output / "active_site_comparison.csv", site_rows)
    _table(output / "all_amide_geometry.csv", all_rows)
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    (output / "README.md").write_text(
        "# Matched original-backbone geometry\n\n"
        "active_site_comparison.csv contains all four physical sites before and after Amber. "
        "Selected rows use the saved one-to-one assignment after positive-phi Tier A/B eligibility; "
        "no new geometry thresholds are applied. Unassigned sites distinguish no eligible amide "
        "from a one-to-one ligand conflict. Nearest-d1 rows are separate diagnostics and may fail the pose gate.\n\n"
        "all_amide_geometry.csv retains every measured site × ligand × amide row (48 per snapshot). "
        "Cases without completed MM results have explicit unavailable site placeholders and no fabricated geometry. "
        "They are excluded from measured-site and pose counts. Original failed evidence remains in the cohort. "
        "JSON is authoritative. Force convergence, pose acceptance and native-Thr similarity remain separate; "
        "these geometric proxies do not establish activity or barriers.\n"
    )
    return summary


def main() -> None:
    """Read an explicit cohort root and write a new comparison directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-root", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output if args.output is not None else args.cohort_root / "comparison"
    summary = summarize_cohort(args.cohort_root, output)
    logging.info(
        "Published %s measured cases; %s cases without MM results",
        summary["measured_mm_cases"],
        summary["cases_without_mm_result"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
