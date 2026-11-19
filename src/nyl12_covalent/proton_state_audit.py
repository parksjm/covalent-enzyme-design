"""Summarize cpptraj distance tables for the window-24 proton-state audit."""

from __future__ import annotations

import json
import logging
from collections import Counter
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)
REQUIRED_COLUMNS = (
    "#Frame",
    "N_HG1",
    "OG1_HG1",
    "N_H1",
    "OG1_H1",
    "N_H2",
    "OG1_H2",
    "OG1_C",
    "C_O",
)
HYDROGENS = ("HG1", "H1", "H2")


def classify_hydrogen(n_distance: float, og1_distance: float) -> str:
    """Classify a hydrogen from its distances to N and OG1.

    Args:
        n_distance: Distance from the hydrogen to the N-terminal nitrogen.
        og1_distance: Distance from the hydrogen to Thr233 OG1.

    Returns:
        ``N-bound``, ``OG1-bound``, or ``ambiguous``.
    """
    if n_distance <= 1.25 and og1_distance >= 1.30:
        return "N-bound"
    if og1_distance <= 1.25 and n_distance >= 1.30:
        return "OG1-bound"
    return "ambiguous"


def _read_table(path: Path) -> list[dict[str, float]]:
    """Read whitespace-delimited cpptraj data and validate its header."""
    with path.open(encoding="utf-8") as stream:
        lines = [line.split() for line in stream if line.strip()]
    if not lines or tuple(lines[0]) != REQUIRED_COLUMNS:
        raise ValueError(f"{path}: expected header {' '.join(REQUIRED_COLUMNS)}")
    records: list[dict[str, float]] = []
    for line_number, fields in enumerate(lines[1:], start=2):
        if len(fields) != len(REQUIRED_COLUMNS):
            raise ValueError(f"{path}:{line_number}: expected 9 fields")
        try:
            records.append(dict(zip(REQUIRED_COLUMNS, map(float, fields), strict=True)))
        except ValueError as error:
            raise ValueError(
                f"{path}:{line_number}: nonnumeric distance row"
            ) from error
    if not records:
        raise ValueError(f"{path}: no distance records")
    return records


def summarize_stage(path: Path) -> dict[str, Any]:
    """Summarize one stage table's frames, hydrogen states, and distances.

    Args:
        path: Whitespace-delimited cpptraj distance table.

    Returns:
        Per-stage state counts and distance statistics.
    """
    records = _read_table(path)
    state_counts: dict[str, dict[str, int]] = {
        hydrogen: Counter() for hydrogen in HYDROGENS
    }
    intended_count = 0
    og1_c_count = 0
    for row in records:
        states = {
            hydrogen: classify_hydrogen(row[f"N_{hydrogen}"], row[f"OG1_{hydrogen}"])
            for hydrogen in HYDROGENS
        }
        for hydrogen, state in states.items():
            state_counts[hydrogen][state] += 1
        intended_count += all(state == "N-bound" for state in states.values())
        og1_c_count += row["OG1_C"] <= 2.5
    distance_stats = {
        column: {
            "min": min(row[column] for row in records),
            "mean": sum(row[column] for row in records) / len(records),
            "max": max(row[column] for row in records),
        }
        for column in REQUIRED_COLUMNS[1:]
    }
    return {
        "stage": path.stem,
        "frame_count": len(records),
        "state_counts": {key: dict(value) for key, value in state_counts.items()},
        "intended_state_frames": intended_count,
        "OG1_C_le_2_5_frames": og1_c_count,
        "distance_stats": distance_stats,
    }


def summarize_tables(paths: list[Path]) -> dict[str, Any]:
    """Summarize multiple stage tables in the given order.

    Args:
        paths: Ordered cpptraj distance tables.

    Returns:
        Mapping containing the ordered stage summaries.
    """
    if not paths:
        raise ValueError("at least one stage table is required")
    stages = [summarize_stage(path) for path in paths]
    return {"stages": stages}


def write_summary(paths: list[Path], output: Path) -> None:
    """Write stage summaries as formatted JSON.

    Args:
        paths: Ordered cpptraj distance tables.
        output: Destination JSON path.
    """
    result = summarize_tables(paths)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    LOGGER.info("Wrote proton-state audit summary to %s", output)
