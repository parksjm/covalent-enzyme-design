"""Rank observed relay paths without changing their chemistry classes."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from nyl12_covalent.modular_es_snapshot_selection import classify_geometry


def rank(row: dict[str, Any]) -> tuple:
    """Return deterministic strict-path, tier and geometric ranking."""
    return (
        not row["strict"],
        {"A": 0, "B": 1, None: 2}[row["tier"]],
        row.get("network_priority", 3),
        row["violation"],
        row["water_count"] + row["internal_count"],
        row["case"],
        row["frame"],
        tuple(row["atom_indices1"]),
    )


def read_candidates(directories: list[Path], config: dict) -> list[dict[str, Any]]:
    """Read best path per distinct frame, allowing refined records to supersede."""
    cases = {case["name"]: case for case in config["cases"]}
    frames: dict[tuple[str, int], dict] = {}
    found = False
    for directory in directories:
        local: dict[tuple[str, int], dict] = {}
        for path in sorted(directory.rglob("*_paths.jsonl")):
            found = True
            with path.open() as handle:
                for line in handle:
                    row = json.loads(line)
                    row["tier"] = classify_geometry(row)
                    phi = row["projected_phi"]
                    row["phi_orientation"] = (
                        "positive" if phi > 0 else "negative" if phi < 0 else "zero"
                    )
                    _annotate_row(row, cases[row["case"]])
                    key = row["case"], row["frame"]
                    if key not in local or rank(row) < rank(local[key]):
                        local[key] = row
        frames.update(local)
    if not found:
        raise ValueError("No *_paths.jsonl screen records found")
    return sorted(frames.values(), key=rank)


def diverse(rows: list[dict], limit: int, separation_ps: float = 1000) -> list[dict]:
    """Prefer independent cases, requiring temporal separation within a case."""
    selected: list[dict] = []
    pending = list(rows)
    while pending and len(selected) < limit:
        used = {r["case"] for r in selected}
        index = next((i for i, row in enumerate(pending) if row["case"] not in used), 0)
        row = pending.pop(index)
        if all(
            row["case"] != prior["case"]
            or abs(row["time_ps"] - prior["time_ps"]) >= separation_ps
            for prior in selected
        ):
            selected.append(row)
    return selected


def shortlist(rows: list[dict], limit: int = 20) -> list[dict]:
    """Bound native evaluation while retaining representative eligible frames."""
    qualified = [r for r in rows if r["strict"] and r["tier"] in {"A", "B"}]
    selected = diverse(qualified, limit)
    keys = {(r["case"], r["frame"]) for r in selected}
    for row in rows:
        key = row["case"], row["frame"]
        if len(selected) == limit:
            break
        if key not in keys:
            selected.append(row)
            keys.add(key)
    return selected


def refinement_frames(rows: list[dict], totals: dict[str, int]) -> list[dict]:
    """Request the ±4 saved-frame neighborhood without exceeding source bounds."""
    return [
        {"case": case, "frame": frame}
        for case, frame in sorted(
            {
                (row["case"], frame)
                for row in rows
                for frame in range(
                    max(1, row["frame"] - 4),
                    min(totals[row["case"]], row["frame"] + 4) + 1,
                )
            }
        )
    ]


def final_selection(rows: list[dict]) -> list[dict]:
    """Select at most two primary and two reserves; keep native failures explicit."""
    eligible = [r for r in rows if r["strict"] and r["tier"] in {"A", "B"}]
    eligible.sort(key=rank)
    selected = diverse(eligible, 4)
    return [
        {**row, "selection_role": "primary" if i < 2 else "reserve"}
        for i, row in enumerate(selected)
    ]


def annotate_context(rows: list[dict], config: dict) -> list[dict]:
    """Label route completeness from explicit starting donor state and mapping."""
    cases = {case["name"]: case for case in config["cases"]}
    for row in rows:
        _annotate_row(row, cases[row["case"]])
    return sorted(rows, key=rank)


def _annotate_row(row: dict, case: dict) -> None:
    """Annotate one route before any destructive per-frame reduction."""
    start = next(role for role in case["roles"] if role["label"] == row["start"])
    og1 = case["geometry_atoms"]["og1"]
    nalpha = case["geometry_atoms"]["nalpha"]

    def identity(descriptor: dict) -> tuple:
        return descriptor["residue_index1"], descriptor["atom"]

    if identity(start) == identity(og1):
        priority, context = 0, "OG1_start_connected_path"
    elif (
        start["state"] in {"NH3+", "COOH", "ASH", "protonated"}
        and not start["conditional"]
    ):
        priority, context = 1, "preprotonated_donor_path"
    else:
        priority, context = 2, "conditional_activation_hypothesis"
    row["network_priority"] = priority
    row["network_context"] = context
    row["nalpha_NH2_start"] = (
        identity(start) == identity(nalpha) and start["state"] == "NH2"
    )
