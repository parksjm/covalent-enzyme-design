"""Parse umbrella coordinates and quantify neighboring-window overlap."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


def read_dumpave(path: Path) -> np.ndarray:
    """Read the sampled coordinate column from an Amber DUMPAVE file.

    Args:
        path: Amber NMR restraint coordinate output.

    Returns:
        One-dimensional sampled coordinate array in angstroms.
    """
    values: list[float] = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 2:
            raise ValueError(f"Malformed DUMPAVE line {line_number} in {path}")
        try:
            int(fields[0])
            values.append(float(fields[1]))
        except ValueError as error:
            raise ValueError(
                f"Non-numeric DUMPAVE line {line_number} in {path}"
            ) from error
    result = np.asarray(values, dtype=float)
    if result.size == 0 or not np.isfinite(result).all():
        raise ValueError(f"No finite umbrella samples in {path}")
    return result


def campaign_samples(
    campaign: Path,
    stage: str = "03_production_2ps",
    discard_fraction: float = 0.0,
    allow_incomplete: bool = False,
) -> list[dict[str, Any]]:
    """Load ordered umbrella samples from a campaign.

    Args:
        campaign: Umbrella campaign containing ``inputs.json`` and windows.
        stage: Explicit stage directory or glob pattern to analyze.
        discard_fraction: Fraction removed from the beginning of each stream.
        allow_incomplete: Permit missing or empty windows for live diagnostics.

    Returns:
        Ordered window records with targets, source paths, and sample arrays.
    """
    if not 0.0 <= discard_fraction < 1.0:
        raise ValueError("discard_fraction must be in [0, 1)")
    inputs = json.loads((campaign / "inputs.json").read_text())
    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for window in inputs["windows"]:
        window_directory = campaign / window["label"]
        stages = (
            sorted(window_directory.glob(stage))
            if any(character in stage for character in "*?[")
            else [window_directory / stage]
        )
        paths: list[Path] = []
        arrays: list[np.ndarray] = []
        for stage_directory in stages:
            process = stage_directory / "process.json"
            path = stage_directory / "umbrella_coordinate.dat"
            if not process.is_file() or not path.is_file():
                continue
            if json.loads(process.read_text()).get("exit_status") != 0:
                continue
            try:
                arrays.append(read_dumpave(path))
            except ValueError:
                continue
            paths.append(path)
        if not arrays:
            if not allow_incomplete:
                missing.append(window["label"])
            continue
        samples = np.concatenate(arrays)
        start = int(samples.size * discard_fraction)
        samples = samples[start:]
        if samples.size == 0:
            if not allow_incomplete:
                missing.append(window["label"])
            continue
        rows.append(
            {
                "window": window["label"],
                "target_A": float(window["target_A"]),
                "source": ";".join(str(path.resolve()) for path in paths),
                "samples_A": samples,
            }
        )
    if missing:
        raise ValueError(f"Incomplete umbrella stage {stage}: {', '.join(missing)}")
    if len(rows) < 2:
        raise ValueError("At least two sampled windows are required")
    return rows


def overlap_summary(
    windows: list[dict[str, Any]], bins: int = 80
) -> tuple[list[dict[str, Any]], list[float]]:
    """Summarize windows and histogram-intersection overlap of neighbors.

    Args:
        windows: Output from :func:`campaign_samples`.
        bins: Number of shared histogram bins.

    Returns:
        Per-window summary rows and common bin edges.
    """
    if len(windows) < 2 or bins < 2:
        raise ValueError("At least two windows and two bins are required")
    combined = np.concatenate([row["samples_A"] for row in windows])
    lower, upper = float(combined.min()), float(combined.max())
    if lower == upper:
        raise ValueError("Umbrella samples have zero coordinate range")
    edges = np.linspace(lower, upper, bins + 1)
    probabilities = []
    summaries: list[dict[str, Any]] = []
    for row in windows:
        counts, _ = np.histogram(row["samples_A"], bins=edges)
        probability = counts / counts.sum()
        probabilities.append(probability)
        summaries.append(
            {
                "window": row["window"],
                "target_A": row["target_A"],
                "sample_count": int(row["samples_A"].size),
                "mean_A": float(np.mean(row["samples_A"])),
                "standard_deviation_A": float(np.std(row["samples_A"], ddof=1))
                if row["samples_A"].size > 1
                else 0.0,
                "minimum_A": float(np.min(row["samples_A"])),
                "maximum_A": float(np.max(row["samples_A"])),
                "neighbor_overlap": None,
                "source": row["source"],
            }
        )
    for index in range(len(summaries) - 1):
        summaries[index]["neighbor_overlap"] = float(
            np.minimum(probabilities[index], probabilities[index + 1]).sum()
        )
    return summaries, edges.tolist()


def write_overlap_data(
    summaries: list[dict[str, Any]], edges: list[float], output: Path
) -> tuple[Path, Path]:
    """Write machine-readable overlap summaries.

    Args:
        summaries: Output rows from :func:`overlap_summary`.
        edges: Common histogram bin edges.
        output: Existing output directory.

    Returns:
        JSON and CSV paths.
    """
    json_path = output / "overlap_summary.json"
    csv_path = output / "overlap_summary.csv"
    json_path.write_text(
        json.dumps({"bin_edges_A": edges, "windows": summaries}, indent=2) + "\n"
    )
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    return json_path, csv_path
