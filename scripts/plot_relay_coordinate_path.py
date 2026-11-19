#!/usr/bin/env python3
"""Plot completed coordinate-driving endpoints; no NEB or minimizer trajectories."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import shutil
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
DEFAULT = (
    ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation"
)


def digest(path: Path) -> str:
    """Return SHA-256 of an input or source artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_complete(directory: Path) -> tuple[dict[str, Any], dict[str, str]]:
    """Read the final complete accepted row of an explicitly completed stage.

    Args:
        directory: Stage containing process.json and accepted_iterations.jsonl.

    Returns:
        Final accepted record and input-file SHA-256 identities.

    Raises:
        ValueError: If stage is incomplete or no complete accepted row exists.
    """
    process, ledger = (
        directory / "process.json",
        directory / "accepted_iterations.jsonl",
    )
    if json.loads(process.read_text()).get("completed") is not True:
        raise ValueError(f"Stage is not complete: {directory}")
    final = None
    with ledger.open() as stream:
        for line in stream:
            if line.endswith("\n"):
                final = json.loads(line)
    if final is None:
        raise ValueError(f"No complete accepted row: {ledger}")
    return final, {str(p): digest(p) for p in (process, ledger)}


def render(campaign: Path, output: Path) -> dict[str, Any]:
    """Write a forward E–S-to-acyl plot, inspection table and authoritative JSON.

    Args:
        campaign: Completed water-relay campaign directory.
        output: Output directory for PNG, PDF, CSV, JSON and source snapshot.

    Returns:
        Numeric summary and input/source identities retained in the JSON output.
    """
    stages = [("ES", campaign / "03_reverse_path/unbiased_ES_release")]
    stages += [
        (f"window_{i:03d}", campaign / f"03_reverse_path/window_{i:03d}")
        for i in range(36, 0, -1)
    ]
    stages.append(("acyl", campaign / "02_acyl/02_unbiased_release"))
    records, hashes = [], {}
    for name, directory in stages:
        row, identities = last_complete(directory)
        records.append((name, row))
        hashes.update(identities)
    baseline = records[0][1]["physical_energy_kcal_mol"]
    table = []
    for point, (name, row) in enumerate(records):
        physical, bias, total = (
            float(row[k])
            for k in (
                "physical_energy_kcal_mol",
                "bias_energy_kcal_mol",
                "total_energy_kcal_mol",
            )
        )
        if not np.isfinite([physical, bias, total]).all() or not np.isclose(
            physical + bias, total, rtol=0, atol=1e-6
        ):
            raise ValueError(f"Invalid/inconsistent energies: {name}")
        table.append(
            {
                "point": point,
                "stage": name,
                "final_iteration": row["iteration"],
                "physical_energy_kcal_mol": physical,
                "bias_energy_kcal_mol": bias,
                "total_energy_kcal_mol": total,
                "relative_physical_kcal_mol": physical - baseline,
                "relative_total_kcal_mol": total - baseline,
            }
        )
    if any(abs(table[i]["bias_energy_kcal_mol"]) > 1e-10 for i in (0, -1)):
        raise ValueError("E–S and acyl endpoints must have zero bias")
    output.mkdir(parents=True, exist_ok=True)
    prefix = output / "coordinate_driving_pes"
    with prefix.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    peak = max(table, key=lambda r: r["relative_physical_kcal_mol"])
    summary = {
        "physical_max_kcal_mol": peak["relative_physical_kcal_mol"],
        "physical_max_stage": peak["stage"],
        "acyl_relative_kcal_mol": table[-1]["relative_physical_kcal_mol"],
    }
    report = {
        "definition": "Final completed accepted state per driving window; unbiased endpoints",
        "direction": "E–S → acyl-enzyme; windows036→001",
        "reference_physical_energy_kcal_mol": baseline,
        "interpretation": "Coordinate-driving PES;25↔24structural jump;not continuous MEP,verified TS,PMF or NEB",
        "source_sha256": {str(Path(__file__).resolve()): digest(Path(__file__))},
        "input_sha256": hashes,
        "summary": summary,
        "points": table,
    }
    prefix.with_suffix(".json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    shutil.copy2(__file__, output / Path(__file__).name)
    x = np.arange(len(table))
    jump = next(i for i, row in enumerate(table) if row["stage"] == "window_025")
    fig, ax = plt.subplots(figsize=(10.2, 5.7))
    for key, color, label, width in (
        ("relative_physical_kcal_mol", "#1764ad", "Physical energy", 2.1),
        ("relative_total_kcal_mol", "#df8526", "Physical + driving bias", 1.1),
    ):
        y = np.array([row[key] for row in table])
        style = {"color": color, "lw": width, "marker": "o", "ms": 3}
        ax.plot(x[: jump + 1], y[: jump + 1], label=label, **style)
        ax.plot(x[jump + 1 :], y[jump + 1 :], color=color, lw=width, marker="o", ms=3)
        ax.plot(x[jump : jump + 2], y[jump : jump + 2], color=color, lw=width, ls="--")
    ax.axvspan(jump, jump + 1, color="0.5", alpha=0.12)
    ticks = [0, 1, 7, 12, 13, 17, 22, 27, 32, 36, 37]
    ax.set_xticks(
        ticks, ["E–S", "36", "30", "25", "24", "20", "15", "10", "5", "1", "Acyl"]
    )
    ax.set(
        xlabel="Reaction-path point (E–S → acyl-enzyme)",
        ylabel="Relative energy (kcal/mol)",
        title="ASP272/ASH274—coordinate driving",
    )
    for point, label, offset in [
        (0, "E–S: 0.00", 6),
        (37, f"Acyl: {summary['acyl_relative_kcal_mol']:.2f}", -6),
    ]:
        ax.annotate(
            label,
            (point, table[point]["relative_physical_kcal_mol"]),
            xytext=(offset, 12),
            textcoords="offset points",
            fontsize=10,
            ha="left" if point == 0 else "right",
        )
    ax.annotate(
        f"{summary['physical_max_kcal_mol']:.2f} ({peak['stage'].replace('_', ' ')})",
        (peak["point"], peak["relative_physical_kcal_mol"]),
        xytext=(8, 10),
        textcoords="offset points",
        fontsize=10,
        color="#1764ad",
    )
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, loc="upper left")
    ax.margins(x=0.025, y=0.16)
    fig.subplots_adjust(bottom=0.22, top=0.9, left=0.1, right=0.98)
    fig.text(
        0.1,
        0.065,
        "Dashed segment/shading: window 25↔24 structural jump; not a continuous MEP or verified TS.\n"
        "Final accepted state per window; zero-bias endpoints. Reference: final E–S energy. No NEB data.",
        fontsize=9,
    )
    fig.savefig(prefix.with_suffix(".png"), dpi=220)
    fig.savefig(prefix.with_suffix(".pdf"))
    plt.close(fig)
    LOGGER.info("Wrote %s.png; %s", prefix, summary)
    return report


def main() -> None:
    """Parse paths and render the completed-path figure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=DEFAULT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    render(
        campaign,
        args.output.resolve() if args.output else campaign / "03_reverse_path/plots",
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
