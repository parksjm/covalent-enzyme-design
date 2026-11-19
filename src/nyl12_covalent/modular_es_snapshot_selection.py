"""Select temporally representative modular-PA6 E-S snapshots."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class GeometryThresholds:
    """Operational near-attack geometry limits."""

    d1_max: float
    oxyanion_max: float
    attack_min: float
    attack_max: float
    abs_phi_min: float
    abs_phi_max: float


TIERS = {
    "A": GeometryThresholds(3.2, 3.8, 90.0, 120.0, 70.0, 110.0),
    "B": GeometryThresholds(3.5, 4.0, 85.0, 125.0, 60.0, 125.0),
}
METRIC_COLUMNS = (
    "d1",
    "d2",
    "d3",
    "attack",
    "ligand_rmsd",
    "projected_phi",
)


def _finite(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid numeric value for {label}: {value!r}") from error
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite numeric value for {label}")
    return number


def classify_geometry(row: dict[str, float]) -> str | None:
    """Return the strictest Tier A/B class passed by one geometry row."""
    for name, limits in TIERS.items():
        if (
            row["d1"] <= limits.d1_max
            and row["d2"] <= limits.oxyanion_max
            and row["d3"] <= limits.oxyanion_max
            and limits.attack_min <= row["attack"] <= limits.attack_max
            and limits.abs_phi_min <= abs(row["projected_phi"]) <= limits.abs_phi_max
        ):
            return name
    return None


def read_geometry_timeseries(path: Path) -> list[dict[str, float]]:
    """Read and validate the modular PA6 geometry-timeseries schema."""
    required = {"frame", "time_ns", *METRIC_COLUMNS}
    rows: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Geometry CSV is missing columns: {sorted(required)}")
        for ordinal, source in enumerate(reader, start=1):
            row = {key: _finite(source[key], key) for key in required}
            if not row["frame"].is_integer() or int(row["frame"]) != ordinal:
                raise ValueError("Frame indices must be sequential one-based integers")
            if row["time_ns"] < 0 or any(row[key] <= 0 for key in ("d1", "d2", "d3")):
                raise ValueError("Times must be nonnegative and distances positive")
            rows.append(row)
    if not rows:
        raise ValueError(f"No geometry rows in {path}")
    if any(right["time_ns"] <= left["time_ns"] for left, right in zip(rows, rows[1:])):
        raise ValueError("Geometry times must increase strictly")
    return rows


def select_modular_es_frames(
    csv_path: Path,
    windows_ns: list[tuple[float, float]],
) -> list[dict[str, Any]]:
    """Select one strictest-tier frame nearest each time-window midpoint.

    Tier A is preferred to Tier B, but geometry does not rank frames within the
    selected tier. Ties resolve to the earlier frame. Windows are half-open and
    must be nonoverlapping.
    """
    if not windows_ns:
        raise ValueError("At least one selection window is required")
    windows = [
        (_finite(a, "window start"), _finite(b, "window end")) for a, b in windows_ns
    ]
    if any(start < 0 or start >= end for start, end in windows):
        raise ValueError("Windows require 0 <= start < end")
    ordered = sorted(windows)
    if any(right[0] < left[1] for left, right in zip(ordered, ordered[1:])):
        raise ValueError("Selection windows overlap")

    rows = read_geometry_timeseries(csv_path)
    selected: list[dict[str, Any]] = []
    for start, end in windows:
        candidates = [row for row in rows if start <= row["time_ns"] < end]
        classified = [(row, classify_geometry(row)) for row in candidates]
        counts = {name: sum(tier == name for _, tier in classified) for name in TIERS}
        tier = next((name for name in TIERS if counts[name]), None)
        if tier is None:
            raise ValueError(f"No Tier A/B frame in [{start}, {end}) ns")
        eligible = [row for row, observed in classified if observed == tier]
        midpoint = (start + end) / 2.0
        chosen = min(
            eligible, key=lambda row: (abs(row["time_ns"] - midpoint), row["frame"])
        )
        selected.append(
            {
                **chosen,
                "frame": int(chosen["frame"]),
                "tier": tier,
                "window_start_ns": start,
                "window_end_ns": end,
                "window_total_count": len(candidates),
                "window_tier_counts": counts,
            }
        )
    return selected
