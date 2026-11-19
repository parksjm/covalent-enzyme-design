"""Analyze matched embedded-QM energies along an enzyme reaction path."""

from __future__ import annotations

import csv
import json
import logging
import math
from pathlib import Path
import re
from typing import Any

LOGGER = logging.getLogger(__name__)
HARTREE_TO_KCAL_MOL = 627.5094740631
DFTB_PATTERN = re.compile(r"DFTBESCF=\s*([-+0-9.EeDd]+)")
EXTERN_PATTERN = re.compile(r"EXTERNESCF=\s*([-+0-9.EeDd]+)")
EPTOT_PATTERN = re.compile(r"EPtot\s*=\s*([-+0-9.EeDd]+)")
EAMBER_PATTERN = re.compile(r"EAMBER \(non-restraint\)\s*=\s*([-+0-9.EeDd]+)")
XTB_PATTERN = re.compile(r"\| TOTAL ENERGY\s+([-+0-9.EeDd]+)\s+Eh")


def _last(path: Path, pattern: re.Pattern[str], termination: str) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    if termination not in text:
        raise ValueError(f"Missing termination marker in {path}")
    matches = pattern.findall(text)
    if not matches:
        raise ValueError(f"Missing energy field in {path}")
    return float(matches[-1].replace("D", "E").replace("d", "e"))


def read_window(window: Path) -> dict[str, float]:
    """Read the three embedded-QM energies and Sander hybrid potentials.

    Args:
        window: One completed benchmark window directory.

    Returns:
        Absolute energies in kcal/mol.
    """
    dftb = window / "dftb3_mm" / "single_point.mdout"
    dft = window / "wb97xd_mm" / "single_point.mdout"
    xtb = window / "gfn2_xtb" / "gfn2.out"
    return {
        "dftb3_qm_kcal_mol": _last(dftb, DFTB_PATTERN, "DFTBESCF="),
        "dftb3_hybrid_kcal_mol": _last(
            dftb, EPTOT_PATTERN, "Energy calculation finished"
        ),
        "gfn2_xtb_qm_kcal_mol": _last(xtb, XTB_PATTERN, "normal termination of xtb")
        * HARTREE_TO_KCAL_MOL,
        "wb97xd_qm_kcal_mol": _last(dft, EXTERN_PATTERN, "EXTERNESCF="),
        "wb97xd_hybrid_kcal_mol": _last(
            dft, EPTOT_PATTERN, "Energy calculation finished"
        ),
    }


def analyze(input_root: Path, output: Path) -> None:
    """Write matched relative-energy tables, metrics, and a profile plot.

    Args:
        input_root: Completed benchmark input/output root.
        output: New analysis-result directory.
    """
    manifest = json.loads((input_root / "benchmark_manifest.json").read_text())
    rows: list[dict[str, Any]] = []
    for record in manifest["records"]:
        row = {"window": record["window"]}
        row.update(record["acceptance"]["distances_A"])
        row.update(read_window(input_root / "windows" / record["window"]))
        row["mm_remainder_kcal_mol"] = (
            row["dftb3_hybrid_kcal_mol"] - row["dftb3_qm_kcal_mol"]
        )
        rows.append(row)
    methods = ("dftb3_qm", "gfn2_xtb_qm", "wb97xd_qm", "mm_remainder")
    for method in methods:
        origin = rows[0][f"{method}_kcal_mol"]
        for row in rows:
            row[f"relative_{method}_kcal_mol"] = row[f"{method}_kcal_mol"] - origin
    for row in rows:
        row["relative_gfn2_xtb_hybrid_kcal_mol"] = (
            row["relative_gfn2_xtb_qm_kcal_mol"] + row["relative_mm_remainder_kcal_mol"]
        )
        row["relative_dftb3_hybrid_kcal_mol"] = (
            row["dftb3_hybrid_kcal_mol"] - rows[0]["dftb3_hybrid_kcal_mol"]
        )
        row["relative_wb97xd_hybrid_kcal_mol"] = (
            row["wb97xd_hybrid_kcal_mol"] - rows[0]["wb97xd_hybrid_kcal_mol"]
        )
    metrics: dict[str, Any] = {}
    for approximate in ("dftb3", "gfn2_xtb"):
        qm_errors = [
            row[f"relative_{approximate}_qm_kcal_mol"]
            - row["relative_wb97xd_qm_kcal_mol"]
            for row in rows
        ]
        hybrid_errors = [
            row[f"relative_{approximate}_hybrid_kcal_mol"]
            - row["relative_wb97xd_hybrid_kcal_mol"]
            for row in rows
        ]
        metrics[approximate] = {
            "qm_endpoint_error_kcal_mol": qm_errors[-1],
            "qm_rmse_kcal_mol": math.sqrt(
                sum(value * value for value in qm_errors) / len(qm_errors)
            ),
            "hybrid_endpoint_error_kcal_mol": hybrid_errors[-1],
            "hybrid_rmse_kcal_mol": math.sqrt(
                sum(value * value for value in hybrid_errors) / len(hybrid_errors)
            ),
        }
    output.mkdir(parents=True, exist_ok=False)
    (output / "matched_energies.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    with (output / "matched_energies.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "method_errors.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )
    _plot(rows, output / "embedded_method_comparison.png")
    LOGGER.info("Analyzed %d matched enzyme-path geometries in %s", len(rows), output)


def correct_periodic_profile(
    matched_energies: Path, coordinate_drive: Path, output: Path
) -> None:
    """Apply embedded DFT-minus-DFTB3 differences to a periodic DFTB3 path.

    The absolute finite-cutoff profiles are not compared across windows because
    Amber's 20-A point-charge membership changes with geometry. The per-window
    DFTB3-minus-DFT double difference is instead applied to the original
    periodic-Ewald DFTB3 QM and non-restraint hybrid profiles.

    Args:
        matched_energies: JSON table from :func:`analyze`.
        coordinate_drive: Original periodic direct-proton-transfer path.
        output: New result directory.
    """
    matched = json.loads(matched_energies.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for source in matched:
        window = source["window"]
        mdout = coordinate_drive / window / "minimized.mdout"
        periodic_qm = _last(mdout, DFTB_PATTERN, "DFTBESCF=")
        periodic_hybrid = _last(mdout, EAMBER_PATTERN, "EAMBER (non-restraint)")
        dftb_error = (
            source["relative_dftb3_qm_kcal_mol"] - source["relative_wb97xd_qm_kcal_mol"]
        )
        rows.append(
            {
                "window": window,
                "proton_difference_A": source["proton_difference_A"],
                "finite_embedding_dftb3_minus_dft_kcal_mol": dftb_error,
                "periodic_dftb3_qm_kcal_mol": periodic_qm,
                "periodic_dftb3_hybrid_kcal_mol": periodic_hybrid,
            }
        )
    for field in ("periodic_dftb3_qm", "periodic_dftb3_hybrid"):
        origin = rows[0][f"{field}_kcal_mol"]
        for row in rows:
            relative = row[f"{field}_kcal_mol"] - origin
            row[f"relative_{field}_kcal_mol"] = relative
            row[f"dft_corrected_{field.removeprefix('periodic_dftb3_')}_kcal_mol"] = (
                relative - row["finite_embedding_dftb3_minus_dft_kcal_mol"]
            )
    interpretation = {
        "correction": (
            "Original periodic DFTB3 relative energy minus the per-window "
            "finite-embedding DFTB3-minus-DFT double difference"
        ),
        "raw_finite_cutoff_profile_accepted": False,
        "reason": (
            "The 20-A point-charge membership changes between windows, causing "
            "arbitrary electrostatic-potential offsets for the net -2 QM region."
        ),
        "stationary_point_claim": False,
        "rows": rows,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "dft_corrected_periodic_profile.json").write_text(
        json.dumps(interpretation, indent=2) + "\n", encoding="utf-8"
    )
    with (output / "dft_corrected_periodic_profile.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    _plot_correction(rows, output / "dft_corrected_periodic_profile.png")
    LOGGER.info("Corrected %d periodic path windows in %s", len(rows), output)


def _plot_correction(rows: list[dict[str, Any]], output: Path) -> None:
    import matplotlib.pyplot as plt

    coordinate = [row["proton_difference_A"] for row in rows]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    axes[0].plot(
        coordinate,
        [row["finite_embedding_dftb3_minus_dft_kcal_mol"] for row in rows],
        "o-",
        color="#6a3d9a",
    )
    axes[0].set_title("DFTB3 error relative to DFT")
    axes[0].set_ylabel("Double difference (kcal/mol)")
    axes[1].plot(
        coordinate,
        [row["relative_periodic_dftb3_hybrid_kcal_mol"] for row in rows],
        "o-",
        color="#377eb8",
        label="Periodic DFTB3/MM",
    )
    axes[1].plot(
        coordinate,
        [row["dft_corrected_hybrid_kcal_mol"] for row in rows],
        "s--",
        color="#111111",
        label="DFT-corrected estimate",
    )
    axes[1].set_title("Original periodic path")
    axes[1].set_ylabel("Relative non-restraint energy (kcal/mol)")
    axes[1].legend()
    for axis in axes:
        axis.set_xlabel("d(OG1-HG1) - d(N-HG1) (A)")
        axis.axhline(0.0, color="#999999", linewidth=0.8)
        axis.grid(alpha=0.2)
    figure.suptitle("Five-window DFT correction of the direct proton-transfer path")
    figure.savefig(output, dpi=220)
    plt.close(figure)


def _plot(rows: list[dict[str, Any]], output: Path) -> None:
    import matplotlib.pyplot as plt

    coordinate = [row["proton_difference_A"] for row in rows]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    labels = {
        "dftb3": "DFTB3/3OB",
        "gfn2_xtb": "GFN2-xTB",
        "wb97xd": "ωB97X-D/def2-SVP",
    }
    colors = {"dftb3": "#377eb8", "gfn2_xtb": "#e41a1c", "wb97xd": "#111111"}
    for method in labels:
        axes[0].plot(
            coordinate,
            [row[f"relative_{method}_qm_kcal_mol"] for row in rows],
            "o-",
            color=colors[method],
            label=labels[method],
        )
        axes[1].plot(
            coordinate,
            [row[f"relative_{method}_hybrid_kcal_mol"] for row in rows],
            "o-",
            color=colors[method],
            label=labels[method],
        )
    axes[0].set_title("Embedded QM terms")
    axes[1].set_title("Hybrid potential; common MM remainder for xTB")
    for axis in axes:
        axis.set_xlabel("d(OG1–HG1) − d(N–HG1) (Å)")
        axis.set_ylabel("Relative energy (kcal/mol)")
        axis.axhline(0.0, color="#999999", linewidth=0.8)
        axis.grid(alpha=0.2)
        axis.legend()
    figure.suptitle("Nyl12 direct proton-transfer five-window method comparison")
    figure.savefig(output, dpi=220)
    plt.close(figure)
