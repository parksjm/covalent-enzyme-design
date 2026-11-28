#!/usr/bin/env python3
"""Freeze and analyze the existing WT/mutant static-model control panel."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
LOGGER = logging.getLogger(__name__)


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a required nonempty source table.

    Args:
        path: Explicit CSV input.

    Returns:
        Source rows without inferring experimental labels.
    """
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Empty source table: {path}")
    return rows


def prepare(config: Path, output: Path) -> None:
    """Retain actual inputs and executed sources before analysis.

    Args:
        config: Explicit source/cohort specification.
        output: New analysis directory; existing artifacts are never overwritten.
    """
    cfg = json.loads(config.read_text())
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    records = {}

    def retain(source: Path, relative: Path) -> None:
        destination = output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        records[str(relative)] = {
            "original": str(source.absolute()),
            "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
        }

    retain(config, Path("inputs/original_config.json"))
    for label, value in cfg["inputs"].items():
        retain(Path(value), Path("inputs") / f"{label}{Path(value).suffix}")
    for relative in [
        "scripts/calibrate_nyl12_filters.py",
        "src/nyl12_covalent/filter_calibration.py",
        "src/nyl12_covalent/filter_calibration_plots.py",
        "tests/test_filter_calibration.py",
        "tests/test_filter_calibration_plots.py",
        "tests/test_filter_calibration_entrypoint.py",
        "pyproject.toml",
        "uv.lock",
    ]:
        retain(ROOT / relative, Path("source_snapshot") / relative)
    cfg["prepared_utc"] = datetime.now(timezone.utc).isoformat()
    cfg["inputs"] = {
        key: f"inputs/{key}{Path(value).suffix}" for key, value in cfg["inputs"].items()
    }
    prepared = output / "prepared.json"
    prepared.write_text(json.dumps(cfg, indent=2) + "\n")
    records["prepared.json"] = {
        "sha256": hashlib.sha256(prepared.read_bytes()).hexdigest()
    }
    (output / "input_source_manifest.json").write_text(
        json.dumps({"files": records}, indent=2) + "\n"
    )
    LOGGER.info("Frozen calibration inputs and sources: %s", output)


def load_models(path: Path, cohort: dict) -> list[dict]:
    """Select the prespecified register and model indices, without score selection.

    Args:
        path: Frozen original model table.
        cohort: Explicit prefix, variant set, budget and relaxation-source contract.

    Returns:
        Normalized model records retaining original identities and scores.
    """
    selected = []
    for raw in read_csv(path):
        if not raw["sequence_name"].startswith(cohort["prefix"]):
            continue
        index = int(raw["model"].rsplit("_model_", 1)[1])
        if index < 0:
            raise ValueError("Negative model index")
        if index >= cohort["expected_models"]:
            continue  # The frozen historical cohort explicitly uses indices 0--49.
        if raw["score_source"] != cohort["score_source"]:
            raise ValueError(f"Mixed preparation sources in {raw['model']}")
        selected.append(
            {
                **raw,
                "variant": raw["sequence_name"].removeprefix(cohort["prefix"]).upper(),
                "model_index": index,
                **{
                    key: float(raw[key])
                    for key in (
                        "d1",
                        "d2",
                        "d3",
                        "a1",
                        "herschlag_phi",
                        "interface_delta",
                    )
                },
            }
        )
    return selected


def write_table(path: Path, rows: list[dict]) -> None:
    """Write a nonempty diagnostic table with explicit nested JSON cells.

    Args:
        path: Destination table.
        rows: Result rows with stable column names.
    """
    if not rows:
        raise ValueError(f"Refusing an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(
            {
                key: json.dumps(value) if isinstance(value, (dict, list)) else value
                for key, value in row.items()
            }
            for row in rows
        )


def analyze(prepared: Path) -> None:
    """Audit frozen controls, retain all failures, and export a descriptive report.

    Args:
        prepared: Retained input/cohort specification from the preparation step.
    """
    from nyl12_covalent.filter_calibration import analyze_cohort
    from nyl12_covalent.filter_calibration_plots import export_calibration

    out = prepared.parent
    cfg = json.loads(prepared.read_text())
    manifest = json.loads((out / "input_source_manifest.json").read_text())
    for relative, item in manifest["files"].items():
        if hashlib.sha256((out / relative).read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Frozen input/source changed: {relative}")
    activities = read_csv(out / cfg["inputs"]["activity"])
    for row in activities:
        row["relative_activity"] = float(row["relative_activity"])
    reference = read_csv(out / cfg["inputs"]["historical_summary"])
    results = []
    for cohort in cfg["cohorts"]:
        records = load_models(out / cfg["inputs"][cohort["model_input"]], cohort)
        result = analyze_cohort(records, activities, cohort)
        for row in result["per_variant"]:
            prior = [
                r
                for r in reference
                if (
                    r["substrate_register"] == cohort["register"]
                    and r["variant"] == row["variant"]
                )
            ]
            if len(prior) != 1 or any(
                int(prior[0][key]) != row[key] for key in ["n_models", "n_tier_ab"]
            ):
                raise ValueError(f"Historical count disagreement: {row['variant']}")
        results.append(result)
    report = {
        "analyzed_utc": datetime.now(timezone.utc).isoformat(),
        "cohorts": results,
        "scope": cfg["scope"],
        "assay_caveats": cfg["assay_caveats"],
        "thresholds_changed": False,
        "activity_predictor_validated": False,
        "prospective_design_generation_authorized": False,
        "input_source_manifest_sha256": hashlib.sha256(
            (out / "input_source_manifest.json").read_bytes()
        ).hexdigest(),
    }
    (out / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    for result in results:
        for kind in ["per_variant", "per_model"]:
            write_table(out / "tables" / f"{result['id']}_{kind}.csv", result[kind])
    mismatches = sum(len(r["tier_ab_mismatches"]) for r in results)
    acceptance = {
        "input_source_hashes_verified": True,
        "historical_counts_reproduced": True,
        "tier_ab_parity_mismatches": mismatches,
        "calibration_for_available_derived_labels_complete": mismatches == 0,
        "activity_predictor_validated": False,
        "threshold_tuning_performed": False,
        "raw_assay_uncertainty_available": False,
    }
    (out / "scientific_acceptance.json").write_text(
        json.dumps(acceptance, indent=2) + "\n"
    )
    if mismatches:
        raise ValueError("Tier parity failed; inspect report before interpretation")
    export_calibration(report, out / "plots")
    LOGGER.info("Calibration complete: %s", out / "report.json")


def main() -> None:
    """Prepare or analyze the prespecified retrospective control panel."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--config", type=Path, required=True)
    prep.add_argument("--output", type=Path, required=True)
    run = sub.add_parser("analyze")
    run.add_argument("--prepared", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.config, args.output)
    else:
        analyze(args.prepared)


if __name__ == "__main__":
    main()
