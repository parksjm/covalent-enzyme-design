"""Compare matched DFTB3 and GFN2-xTB profiles with reference DFT."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
from collections import defaultdict
from pathlib import Path
import re
from typing import Any

LOGGER = logging.getLogger(__name__)
HARTREE_TO_KCAL_MOL = 627.5094740631
DFTB_PATTERN = re.compile(r"Total SCF energy\s*=\s*([-+0-9.Ee]+)\s+kcal/mol")
XTB_PATTERN = re.compile(r"\| TOTAL ENERGY\s+([-+0-9.Ee]+)\s+Eh")
DFT_PATTERN = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s*([-+0-9.Ee]+)")
METHODS = ("dftb3", "gfn2_xtb", "wb97xd_augccpvtz")
APPROXIMATE_METHODS = METHODS[:2]


def _single_match(path: Path, pattern: re.Pattern[str], termination: str) -> float:
    """Extract the final numeric match and require normal engine termination."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if termination not in text:
        raise RuntimeError(f"Missing normal termination marker: {path}")
    matches = pattern.findall(text)
    if not matches:
        raise ValueError(f"Energy field not found: {path}")
    return float(matches[-1])


def read_energies(calculation: Path) -> dict[str, float]:
    """Read matched method energies from one calculation directory.

    Args:
        calculation: Directory containing the three engine outputs.

    Returns:
        Absolute electronic energies in kcal/mol.
    """
    return {
        "dftb3": _single_match(
            calculation / "dftb3.out",
            DFTB_PATTERN,
            "Calculation Completed",
        ),
        "gfn2_xtb": _single_match(
            calculation / "gfn2.out",
            XTB_PATTERN,
            "normal termination of xtb",
        )
        * HARTREE_TO_KCAL_MOL,
        "wb97xd_augccpvtz": _single_match(
            calculation / "wb97xd_augccpvtz.log",
            DFT_PATTERN,
            "Normal termination of Gaussian 16",
        )
        * HARTREE_TO_KCAL_MOL,
    }


def _profile_metrics(rows: list[dict[str, Any]], method: str) -> dict[str, float]:
    """Calculate profile errors relative to DFT for one approximate method."""
    errors = [row[f"relative_{method}"] - row["relative_wb97xd_augccpvtz"] for row in rows]
    reference_barrier = max(row["relative_wb97xd_augccpvtz"] for row in rows)
    method_barrier = max(row[f"relative_{method}"] for row in rows)
    return {
        "signed_mean_error_kcal_mol": sum(errors) / len(errors),
        "mean_unsigned_error_kcal_mol": sum(abs(error) for error in errors) / len(errors),
        "rmse_kcal_mol": math.sqrt(sum(error * error for error in errors) / len(errors)),
        "maximum_absolute_error_kcal_mol": max(abs(error) for error in errors),
        "endpoint_reaction_energy_error_kcal_mol": errors[-1],
        "fixed_geometry_maximum_error_kcal_mol": method_barrier - reference_barrier,
        "reference_fixed_geometry_maximum_kcal_mol": reference_barrier,
        "method_fixed_geometry_maximum_kcal_mol": method_barrier,
    }


def _plot_profiles(
    grouped: dict[tuple[str, float], list[dict[str, Any]]],
    manifest: dict[str, Any],
    output: Path,
) -> None:
    """Plot all matched profiles using the analysis environment's Matplotlib."""
    import matplotlib.pyplot as plt

    pairs = [system["name"] for system in manifest["systems"]]
    distances = manifest["distances_A"]
    figure, axes = plt.subplots(len(pairs), len(distances), figsize=(13, 10), sharex=False)
    colors = {
        "dftb3": "#377eb8",
        "gfn2_xtb": "#e41a1c",
        "wb97xd_augccpvtz": "#111111",
    }
    labels = {
        "dftb3": "DFTB3/3OB",
        "gfn2_xtb": "GFN2-xTB",
        "wb97xd_augccpvtz": "ωB97X-D/aug-cc-pVTZ",
    }
    for row_index, pair in enumerate(pairs):
        for column_index, distance in enumerate(distances):
            axis = axes[row_index][column_index]
            profile = grouped[(pair, distance)]
            coordinate = [row["donor_H_A"] - row["acceptor_H_A"] for row in profile]
            for method in METHODS:
                axis.plot(
                    coordinate,
                    [row[f"relative_{method}"] for row in profile],
                    marker="o",
                    color=colors[method],
                    label=labels[method],
                )
            axis.axhline(0.0, color="#999999", linewidth=0.7)
            axis.set_title(f"{pair}; R={distance:.1f} Å")
            axis.set_xlabel("d(donor–H) − d(acceptor–H) (Å)")
            axis.set_ylabel("Relative energy (kcal/mol)")
    handles, legend_labels = axes[0][0].get_legend_handles_labels()
    figure.legend(handles, legend_labels, loc="upper center", ncol=3)
    figure.suptitle("Matched fixed-geometry proton-transfer profiles", y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output / "proton_transfer_method_comparison.png", dpi=220)
    plt.close(figure)


def analyze(input_root: Path, output: Path) -> None:
    """Parse all matched outputs and write tables, metrics, and profile plots.

    Args:
        input_root: Completed generated benchmark directory.
        output: New analysis result directory.
    """
    input_root = input_root.resolve(strict=True)
    manifest = json.loads((input_root / "benchmark_manifest.json").read_text())
    output.mkdir(parents=True, exist_ok=False)
    rows: list[dict[str, Any]] = []
    for record in manifest["records"]:
        calculation = input_root / record["directory"]
        row = dict(record)
        row.update(read_energies(calculation))
        rows.append(row)

    grouped: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(row["pair"], row["donor_acceptor_A"])].append(row)

    metrics: dict[str, Any] = {}
    for (pair, distance), profile in grouped.items():
        profile.sort(key=lambda row: row["point"])
        for method in METHODS:
            origin = profile[0][method]
            for row in profile:
                row[f"relative_{method}"] = row[method] - origin
        key = f"{pair}/R{distance:.1f}"
        metrics[key] = {
            method: _profile_metrics(profile, method)
            for method in APPROXIMATE_METHODS
        }

    fieldnames = [
        "id",
        "pair",
        "charge",
        "donor_acceptor_A",
        "point",
        "donor_H_A",
        "acceptor_H_A",
        *METHODS,
        *(f"relative_{method}" for method in METHODS),
    ]
    with (output / "matched_energies.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    (output / "profile_error_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )

    _plot_profiles(grouped, manifest, output)
    LOGGER.info("Analyzed %d matched points; results: %s", len(rows), output)


def main() -> None:
    """Parse explicit input and output directories."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    analyze(args.input_root, args.output)


if __name__ == "__main__":
    main()
