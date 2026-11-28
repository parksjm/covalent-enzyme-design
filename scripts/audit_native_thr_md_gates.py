#!/usr/bin/env python3
"""Save descriptive original Thr-gate failure attribution and frame selections."""

from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.native_thr_md_gate_audit import (  # noqa: E402
    KEYS,
    LIMITS,
    representatives,
    summarize,
    values,
)


def digest(path: Path) -> str:
    """Return file SHA-256."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Audit saved per-site scalar rows; do not read trajectory payloads."""
    base = ROOT / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison"
    p = base / "md_followup/analysis"
    out = p / "native_thr_md_gate_audit_20261124"
    out.mkdir(exist_ok=True)
    report = {
        "original_limits": dict(zip(KEYS, LIMITS.tolist(), strict=True)),
        "interpretation": "Diagnostic sensitivities only, not proposed new thresholds; sites are nested within one replica; native gate is not an activity assay.",
        "datasets": {},
    }
    selections = {}
    flat = []
    for name in ("WT", "broader_mask07"):
        path = p / f"preliminary_{name}/results/sites.csv"
        rows = list(csv.DictReader(path.open()))
        x = values(rows)
        for i, r in enumerate(rows):
            original = all(x[i] <= LIMITS)
            if original != (r["native_pass"] == "True") or (
                original and r["tier"] in ("A", "B")
            ) != (r["joint_pass"] == "True"):
                raise ValueError("Saved gate mismatch")
        if len(rows) != 4000:
            raise ValueError("Expected 4000 site frames")
        by = {}
        for site in ("all", "AE", "BF", "CG", "DH"):
            block = [r for r in rows if site == "all" or r["site"] == site]
            for cohort, subset in (
                ("all", block),
                ("tier_ab", [r for r in block if r["tier"] in ("A", "B")]),
            ):
                result = summarize(subset)
                by[f"{site}:{cohort}"] = result
                flat.append(
                    {
                        "dataset": name,
                        "site": site,
                        "cohort": cohort,
                        "n": result["n"],
                        **{f"{k}_fails": v for k, v in result["fails"].items()},
                        **{
                            f"{k}_fail_percent": v
                            for k, v in result["fail_percent"].items()
                        },
                        **{
                            f"{k}_exclusive": v
                            for k, v in result["exclusive_fails"].items()
                        },
                        **{
                            f"{k}_leave_one_out_passes": v
                            for k, v in result["leave_one_out_passes"].items()
                        },
                    }
                )
        baseline_path = base / "validation/analysis/occupied_sites.csv"
        selected = json.loads((base / "md_followup/selection.json").read_text())[
            "source_selection"
        ][name]
        baseline = [
            r
            for r in csv.DictReader(baseline_path.open())
            if r["path"] == str(Path(selected) / "relaxed.pdb")
            and r["stage"] == "relaxed"
        ]
        if len(baseline) != 4:
            raise ValueError("Missing four original minimized sites")
        report["datasets"][name] = {
            "source_csv": str(path),
            "source_sha256": digest(path),
            "baseline_sha256": digest(baseline_path),
            "initial_minimized": [
                {
                    "site": r["site"],
                    **{key: float(r[key]) for key in (*KEYS, "scaffold_rmsd")},
                    "native_pass": r["native_pass"] == "True",
                }
                for r in baseline
            ],
            "groups": by,
        }
        chosen = representatives(rows)
        selections[name] = {}
        for label, row in chosen.items():
            time = float(row["time_ps"])
            stage = (
                "05_production_10ns.nc" if time <= 10500 else "06_production_90ns.nc"
            )
            offset = 500 if stage.startswith("05") else 10500
            frame = round((time - offset) / 10) - 1
            if (
                frame < 0
                or abs((time - offset) / 10 - round((time - offset) / 10)) > 1e-6
            ):
                raise ValueError("Frame time does not match trajectory spacing")
            selections[name][label] = {
                "rule": "AE Tier A/B median normalized-gate Euclidean medoid"
                if label == "AE_median_medoid"
                else "Tier A/B minimum maximum normalized gate value; ties time then site",
                "time_ps": time,
                "site": row["site"],
                "trajectory_stage": stage,
                "trajectory_path": str(
                    base / "md_followup/md" / name / "replica_01/dynamics" / stage
                ),
                "frame_index_zero_based": frame,
                "metrics": {
                    k: float(row[k])
                    for k in (
                        *KEYS,
                        "scaffold_rmsd",
                        "d1",
                        "d2",
                        "d3",
                        "a_bd",
                        "herschlag_phi",
                    )
                },
                "native_pass": row["native_pass"] == "True",
                "source_csv_sha256": digest(path),
            }
    (out / "scalar_gate_report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    (out / "selection.json").write_text(
        json.dumps(selections, indent=2, allow_nan=False) + "\n"
    )
    with (out / "scalar_gate_report.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(flat[0]))
        w.writeheader()
        w.writerows(flat)


if __name__ == "__main__":
    main()
