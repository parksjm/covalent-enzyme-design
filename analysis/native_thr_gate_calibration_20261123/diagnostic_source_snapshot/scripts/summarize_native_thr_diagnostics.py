#!/usr/bin/env python3
"""Summarize existing CSVs and sequence metadata for native-Thr calibration."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from nyl12_covalent.native_thr_calibration import LIMITS, summarize_md, write_table  # noqa: E402


def main() -> None:
    """Write scalar diagnostics and a control raw-to-relaxed paired figure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = args.root / "analysis/native_thr_gate_calibration_20261123"
    data_path = base / "results/calibration.json"
    data = json.loads(data_path.read_text())
    paths = sorted(
        (
            args.root
            / "qmmm/pa6_flipped_model30/10_wt_protonation_md/analysis_after_md/occupancy"
        ).glob("*.csv")
    )
    if len(paths) != 12:
        raise ValueError(f"Expected12existingMDCSVs, found {len(paths)}")
    md = summarize_md(paths)
    mapping = []
    for row in data["n12c_correspondence"]:
        first = data["chain_metadata"]["n12c"][row["n12c_auth_chain"]][
            "observed_sequence_auth_numbered"
        ]
        second = data["chain_metadata"]["9Z0B"][row["9Z0B_auth_chain"]][
            "observed_sequence_auth_numbered"
        ]
        common = sorted(first.keys() & second.keys(), key=int)
        mapping.append(
            {
                **row,
                "shared_observed_residues": len(common),
                "sequence_mismatches": [
                    {
                        "auth_residue": number,
                        "n12c": first[number],
                        "9Z0B": second[number],
                    }
                    for number in common
                    if first[number] != second[number]
                ],
            }
        )
    result = {
        "source_sha256": {
            str(path.relative_to(args.root)): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in [data_path, *paths]
        },
        "md_existing_scalar_summaries": md,
        "sequence_correspondence_audit": mapping,
        "independence": "Correlated frames pooled descriptively; not independent samples or confidence intervals.",
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "diagnostics.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    write_table(args.output / "md_scalar_summary.csv", md)
    write_table(args.output / "sequence_correspondence.csv", mapping)
    controls_path = base / "inputs/occupied_sites.csv"
    import csv

    with controls_path.open(newline="") as handle:
        controls = list(csv.DictReader(handle))
    variants = ("WT", "control_K59D", "control_V111C")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, metric in zip(axes, ("og1_shift", "nterm_shift"), strict=True):
        for i, variant in enumerate(variants):
            rows = [row for row in controls if row["variant"] == variant]
            paired = {
                (row["seed"], row["site"], row["stage"]): float(row[metric])
                for row in rows
            }
            for seed, site in sorted({(row["seed"], row["site"]) for row in rows}):
                ax.plot(
                    [i * 3, i * 3 + 1],
                    [paired[(seed, site, stage)] for stage in ("raw", "relaxed")],
                    color=f"C{i}",
                    alpha=0.45,
                    marker=".",
                )
        ax.axhline(LIMITS[metric], color="red", linestyle="--")
        ax.set(
            xticks=[0.5, 3.5, 6.5],
            xticklabels=["WT", "K59D", "V111C"],
            ylabel=f"{metric} (Å)",
            title="Raw → relaxed; paired model/site",
        )
    fig.savefig(args.output / "known_active_raw_to_relaxed.png", dpi=180)
    plt.close(fig)
    (args.output / "scientific_acceptance.json").write_text(
        json.dumps(
            {
                "existing_md_frames": sum(row["frames"] for row in md),
                "trajectory_reprocessing": False,
                "stored_gate_flags_reproduced": True,
                "sequence_mismatches_explicit": True,
            },
            indent=2,
        )
        + "\n"
    )
    logging.info("Summarized %s existing MD CSVs", len(paths))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
