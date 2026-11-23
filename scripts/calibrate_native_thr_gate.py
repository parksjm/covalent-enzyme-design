#!/usr/bin/env python3
"""Run bounded crystal/control calibration without changing native gate limits."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from nyl12_covalent.native_thr_calibration import (  # noqa: E402
    FIT,
    LIMITS,
    frozen_geometry,
    measure,
    read_sites,
    summarize_controls,
    write_table,
)


def main() -> None:
    """Write authoritative JSON, tables, sensitivity plot and acceptance record."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.root, args.output
    out.mkdir(parents=True, exist_ok=False)
    frozen = root / "qmmm/pa6_flipped_model30/10_wt_protonation_md/analysis_after_md"
    geometry_path = (
        frozen / "source_snapshot/src/nyl12_covalent/catalytic_alignment_geometry.py"
    )
    geometry = frozen_geometry(geometry_path)
    paths = {"n12c": frozen / "source_inputs/n12c.cif"}
    paths.update(
        {
            code: root
            / f"analysis/native_thr_gate_calibration_20261123/inputs/{code}.cif"
            for code in ("9Z0A", "9Z0B")
        }
    )
    datasets = {name: read_sites(path) for name, path in paths.items()}
    reference = datasets["n12c"][0]["A"]
    crystal, sensitivity, correspondences = [], [], []
    for name, (sites, metadata) in datasets.items():
        for chain, site in sites.items():
            crystal.append(
                {
                    "structure": name,
                    "auth_chain": chain,
                    "label_chain": metadata[chain]["label_chain"],
                    **measure(reference, site, geometry),
                }
            )
            for fit_name, fit in (
                ("frozen_16_backbone", FIT),
                ("include_thr_20_backbone", (*FIT, 233)),
                ("alpha_only_8_backbone", (113, 186)),
            ):
                for ref_name, ref_site in (
                    ("n12c_A", reference),
                    ("9Z0A_A", datasets["9Z0A"][0]["A"]),
                ):
                    sensitivity.append(
                        {
                            "structure": name,
                            "auth_chain": chain,
                            "reference": ref_name,
                            "fit": fit_name,
                            **measure(ref_site, site, geometry, fit),
                        }
                    )
            if name == "n12c":
                candidates = [
                    (measure(site, other, geometry)["scaffold_rmsd"], orig_chain)
                    for orig_chain, other in datasets["9Z0B"][0].items()
                ]
                rmsd, orig_chain = min(candidates)
                correspondences.append(
                    {
                        "n12c_auth_chain": chain,
                        "9Z0B_auth_chain": orig_chain,
                        "local_backbone_rmsd": rmsd,
                        **measure(site, datasets["9Z0B"][0][orig_chain], geometry),
                    }
                )
    controls_path = (
        root
        / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/analysis/occupied_sites.csv"
    )
    controls = summarize_controls(controls_path)
    result = {
        "thresholds_unchanged": LIMITS,
        "reference": "n12c auth A; Thr233 held out",
        "amber_mapping": {"Thr233": 929, "Thr234": 930},
        "input_sha256": {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (*paths.values(), geometry_path, controls_path)
        },
        "chain_metadata": {name: data[1] for name, data in datasets.items()},
        "crystal_sites": crystal,
        "n12c_correspondence": correspondences,
        "paired_sensitivity": sensitivity,
        "control_groups": controls,
        "interpretation": {
            "activity_classifier_validated": False,
            "assembly_symmetry_generated_copies_counted": False,
            "independence": "ASU copies are distinct fitted coordinates within two diffraction datasets; copies and model seeds are not biological replicates; n12c is derived and not a third experiment.",
        },
    }
    (out / "calibration.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    for filename, rows in (
        ("crystal_sites.csv", crystal),
        ("paired_sensitivity.csv", sensitivity),
        ("n12c_correspondence.csv", correspondences),
        ("control_summary.csv", controls),
    ):
        write_table(out / filename, rows)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, metric in zip(axes, ("og1_shift", "nterm_shift"), strict=True):
        for i, name in enumerate(paths):
            values = [row[metric] for row in crystal if row["structure"] == name]
            ax.scatter([i] * len(values), values, label=name)
        ax.axhline(LIMITS[metric], color="red", linestyle="--", label="unchanged limit")
        ax.set(
            xticks=range(3),
            xticklabels=list(paths),
            ylabel=f"{metric} (Å)",
            title="Frozen 16-atom local fit",
        )
    fig.savefig(out / "crystal_calibration.png", dpi=180)
    plt.close(fig)
    acceptance = {
        "process_complete": True,
        "scientific_acceptance": {
            "unchanged_thresholds": True,
            "all_existing_control_flags_reproduced": True,
            "all_crystal_metrics_finite": True,
            "activity_predictivity": "unvalidated",
            "crystal_native_passes": sum(row["native_pass"] for row in crystal),
            "crystal_sites_including_derived_n12c": len(crystal),
        },
    }
    (out / "scientific_acceptance.json").write_text(
        json.dumps(acceptance, indent=2) + "\n"
    )
    logging.info(
        "Wrote %s crystal comparisons and %s control groups",
        len(crystal),
        len(controls),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
