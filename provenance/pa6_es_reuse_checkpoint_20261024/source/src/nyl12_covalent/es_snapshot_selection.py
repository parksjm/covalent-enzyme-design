"""Select temporally representative ES frames using historical Tier B geometry.

Canonical full-stride10 CSV frame 1 is full trajectory frame 1, at 0.05 ns;
its legacy time estimate is zero. Preserve both times and their offset.
"""

import csv
import math
from pathlib import Path
from typing import Any

D1_MAX_ANGSTROM = 3.5
D2_MAX_ANGSTROM = 4.0
D3_MAX_ANGSTROM = 4.0
BURGI_DUNITZ_RANGE_DEGREES = (85.0, 125.0)
ABS_PHI_RANGE_DEGREES = (60.0, 125.0)
METRIC_COLUMNS = (
    "d1_d_OG1_scissile_C",
    "d2_d_Phe113N_scissile_O",
    "d3_d_Asn186ND2_scissile_O",
    "a1_alpha_Burgi_Dunitz",
    "Herschlag_phi",
)


def _finite_float(value: Any, label: str) -> float:
    """Parse a finite numeric CSV value."""
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid numeric value for {label}: {value!r}") from error
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite value for {label}")
    return number


def _eligible(row: dict[str, Any]) -> bool:
    """Apply the explicitly fixed historical Tier B limits."""
    d1, d2, d3, angle, phi = (row[key] for key in METRIC_COLUMNS)
    return (
        d1 <= D1_MAX_ANGSTROM
        and d2 <= D2_MAX_ANGSTROM
        and d3 <= D3_MAX_ANGSTROM
        and BURGI_DUNITZ_RANGE_DEGREES[0] <= angle <= BURGI_DUNITZ_RANGE_DEGREES[1]
        and ABS_PHI_RANGE_DEGREES[0] <= abs(phi) <= ABS_PHI_RANGE_DEGREES[1]
    )


def select_es_frames(
    csv_path: Path,
    model: str,
    replica: str,
    windows_ns: list[tuple[float, float]],
    full_frame_stride: int,
    full_frame_dt_ns: float,
) -> list[dict[str, Any]]:
    """Choose one eligible frame closest to each physical-time window midpoint.

    Windows are nonoverlapping half-open intervals. Geometry eligibility uses
    fixed historical Tier B limits: d1 <= 3.5 A, d2/d3 <= 4 A, Burgi-Dunitz
    85--125 degrees and absolute phi 60--125 degrees. Geometry never ranks
    eligible frames. Ties resolve to the earlier frame.

    Args:
        csv_path: Canonical active-site geometry timeseries CSV.
        model: Exact model identifier.
        replica: Exact replica identifier.
        windows_ns: Physical-time windows in caller-desired output order.
        full_frame_stride: Full trajectory frames per CSV sample; canonical 10.
        full_frame_dt_ns: Time between full trajectory frames; canonical 0.05 ns.

    Returns:
        Selected rows with metrics, explicit ordinal/full frame indices, physical
        and logged times, their offset, and total/eligible counts per window.

    Raises:
        ValueError: Inputs, matching rows, sequential indices or windows are
            invalid, or any requested window contains no eligible frame.
    """
    if (
        isinstance(full_frame_stride, bool)
        or not isinstance(full_frame_stride, int)
        or full_frame_stride < 1
    ):
        raise ValueError("full_frame_stride must be a positive integer")
    interval = _finite_float(full_frame_dt_ns, "full_frame_dt_ns")
    if interval <= 0 or not model or not replica or not windows_ns:
        raise ValueError("Positive frame interval, model, replica and windows required")
    windows = []
    for start, end in windows_ns:
        start = _finite_float(start, "window start")
        end = _finite_float(end, "window end")
        if start < 0 or start >= end:
            raise ValueError("Windows require 0 <= start < end")
        windows.append((start, end))
    ordered = sorted(windows)
    if any(right[0] < left[1] for left, right in zip(ordered, ordered[1:])):
        raise ValueError("Windows overlap")

    rows: list[dict[str, Any]] = []
    required = {"model", "replica", "trajectory_kind", "frame", "time_ns_estimate"}
    required.update(METRIC_COLUMNS)
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("CSV is missing required columns")
        for source in reader:
            if (
                source["model"] != model
                or source["replica"] != replica
                or source["trajectory_kind"] != "full"
            ):
                continue
            frame = _finite_float(source["frame"], "frame")
            if not frame.is_integer() or frame != len(rows) + 1:
                raise ValueError("Matching frame indices must be sequential from 1")
            ordinal = int(frame)
            full_frame = (ordinal - 1) * full_frame_stride + 1
            physical_time = full_frame * interval
            logged_time = _finite_float(source["time_ns_estimate"], "time_ns_estimate")
            row: dict[str, Any] = dict(source)
            row.update({key: _finite_float(source[key], key) for key in METRIC_COLUMNS})
            if any(row[key] <= 0 for key in METRIC_COLUMNS[:3]):
                raise ValueError("Distances must be positive")
            row.update(
                frame=ordinal,
                ordinal_frame=ordinal,
                full_frame=full_frame,
                physical_time_ns=physical_time,
                logged_time_ns=logged_time,
                time_offset_ns=physical_time - logged_time,
            )
            rows.append(row)
    if not rows:
        raise ValueError("No matching full-trajectory rows")

    selected = []
    for start, end in windows:
        candidates = [row for row in rows if start <= row["physical_time_ns"] < end]
        eligible = [row for row in candidates if _eligible(row)]
        if not eligible:
            raise ValueError(f"No eligible ES frame in [{start}, {end}) ns")
        midpoint = (start + end) / 2
        chosen = min(
            eligible,
            key=lambda row: (abs(row["physical_time_ns"] - midpoint), row["frame"]),
        )
        selected.append(
            dict(
                chosen,
                window_start_ns=start,
                window_end_ns=end,
                eligible_count=len(eligible),
                total_count=len(candidates),
            )
        )
    return selected
