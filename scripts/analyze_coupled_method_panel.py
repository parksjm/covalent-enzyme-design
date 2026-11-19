"""Validate engine completion and compare embedded method double differences."""

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def analyze(root: Path) -> None:
    """Write relative method errors, completion checks and a scientific figure."""
    data = json.loads((root / "work/energies.json").read_text())
    manifest = json.loads((root / "work/manifest.json").read_text())
    if [r["state"] for r in data] != [r["window"] for r in manifest["records"]]:
        raise ValueError("Incomplete or reordered benchmark states")
    job = json.loads((root / "job_attempt1/status.json").read_text())
    if job["exit_status"] != 0:
        raise ValueError("Benchmark process did not complete successfully")
    methods = ["dftb3_mm", "pm6_dhplus_mm", "gfn2_xtb"]
    origin = data[0]["energies_kcal_mol"]
    rows = []
    hashes = {}
    for state, source in zip(data, manifest["records"], strict=True):
        directory = root / "work/windows" / state["state"]
        for method in ["dftb3_mm", "pm6_dhplus_mm", "wb97xd_mm"]:
            path = directory / method / "single_point.mdout"
            if "Energy calculation finished" not in path.read_text():
                raise ValueError(f"Incomplete engine output: {path}")
            hashes[str(path.relative_to(root))] = hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
        logs = list((directory / "wb97xd_mm").glob("gau_job.log.*"))
        if not logs or not all(
            "Normal termination of Gaussian" in p.read_text() for p in logs
        ):
            raise ValueError("Gaussian normal termination missing")
        if (
            "normal termination of xtb"
            not in (directory / "gfn2_xtb/gfn2.out").read_text()
        ):
            raise ValueError("xTB normal termination missing")
        embedding = json.loads(
            (directory / "gfn2_xtb/xtb_input_manifest.json").read_text()
        )
        row = {
            "state": state["state"],
            "point_charge_count": embedding["point_charge_count"],
            "OG1_C_A": source["geometry_A"]["OG1_carbonyl_C_A"],
        }
        energies = state["energies_kcal_mol"]
        for method in methods:
            row[method] = (energies[method] - origin[method]) - (
                energies["wb97xd_mm"] - origin["wb97xd_mm"]
            )
        rows.append(row)
    metrics = {}
    for method in methods:
        errors = np.array([r[method] for r in rows])
        metrics[method] = {
            "rmse_kcal_mol": float(np.sqrt(np.mean(errors**2))),
            "max_abs_error_kcal_mol": float(np.max(np.abs(errors))),
            "ti1_vs_es_error_kcal_mol": float(errors[-1]),
            "ti1_vs_coupled_pt_error_kcal_mol": float(errors[-1] - errors[2]),
        }
    output = root / "analysis"
    output.mkdir(exist_ok=False)
    result = {
        "rows": rows,
        "metrics": metrics,
        "source_hashes": hashes,
        "interpretation": "Positive error overestimates relative energy versus DFT. Finite-cutoff raw profiles are not reaction energies. No gradient or stationary-point validation.",
    }
    (output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    with (output / "comparison.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    for method, label in zip(
        methods, ["DFTB3/3OB", "PM6-DH+", "GFN2-xTB"], strict=True
    ):
        ax.plot(range(len(rows)), [r[method] for r in rows], "o-", label=label)
    ax.set_xticks(
        range(len(rows)),
        ["E–S", "Direct PT", "Coupled PT", "Attack", "Late attack", "TI1-like"],
    )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Relative-energy error vs DFT (kcal/mol)")
    ax.set_title("Matched embedded energies; errors referenced to E–S")
    ax.legend()
    fig.savefig(output / "method_errors.png", dpi=180)
    plt.close(fig)
    logging.info("Completed benchmark analysis: %s", metrics)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    logging.basicConfig(level=logging.INFO)
    analyze(parser.parse_args().root)
