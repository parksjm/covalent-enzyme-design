"""Review assisted relay hypotheses separately from zero-intermediate controls."""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

from nyl12_covalent.leaving_group_relay_export import digest
from nyl12_covalent.leaving_group_relay_selection import annotate_context, diverse, rank
from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

CLASSES = (
    "joint_strict_tier_ab",
    "tier_ab_orientation_nearpass",
    "strict_non_tier",
    "non_tier_nearpass",
)


def assisted_class(row: dict) -> str:
    """Keep direct controls and each assisted acceptance class explicit."""
    if row["water_count"] + row["internal_count"] < 1:
        return "direct_control"
    if row["tier"] in {"A", "B"}:
        return CLASSES[0] if row["strict"] else CLASSES[1]
    return CLASSES[2] if row["strict"] else CLASSES[3]


def assisted_rank(row: dict) -> tuple:
    """Rank acceptance classes first and failed orientation violation within class."""
    return (CLASSES.index(row["source_class"]), row["violation"], rank(row))


def read_assisted(
    directories: list[Path],
    config: dict,
    retained: Path,
) -> tuple[list[dict], dict, dict[str, int]]:
    """Retain every input row and select best assisted route per frame.

    Args:
        directories: Complete screen output trees in coarse-to-refined order.
        config: Explicit case configuration used to annotate donor context.
        retained: New JSONL destination retaining all annotated raw records.

    Returns:
        Ranked best assisted frames, counts/provenance and case frame totals.
    """
    frames, totals, summaries, source_hashes = {}, {}, {}, {}
    counts: Counter = Counter()
    with retained.open("x") as output:
        for directory in directories:
            local = {}
            for path in sorted(directory.rglob("*_paths.jsonl")):
                summary_path = path.parent / "summary.json"
                if not summary_path.is_file():
                    raise ValueError(f"Incomplete screen lacks summary: {path.parent}")
                if str(summary_path) not in summaries:
                    summary = json.loads(summary_path.read_text())
                    summaries[str(summary_path)] = digest(summary_path)
                    for case in summary["cases"]:
                        totals[case["case"]] = case["total_frames"]
                        source_hashes.update(case["input_sha256"])
                with path.open() as handle:
                    for number, line in enumerate(handle, 1):
                        row = json.loads(line)
                        if row["case"] not in totals:
                            raise ValueError(
                                "Row case missing from completed summaries"
                            )
                        row["tier"] = classify_geometry(row)
                        row["phi_orientation"] = (
                            "positive"
                            if row["projected_phi"] > 0
                            else "negative"
                            if row["projected_phi"] < 0
                            else "zero"
                        )
                        annotate_context([row], config)
                        row["source_class"] = assisted_class(row)
                        row["source_jsonl"] = str(path.resolve())
                        row["source_line"] = number
                        output.write(json.dumps(row, allow_nan=False) + "\n")
                        counts[row["source_class"]] += 1
                        if row["source_class"] == "direct_control":
                            counts["direct_strict_tier_ab_paths"] += int(
                                row["strict"] and row["tier"] in {"A", "B"}
                            )
                            continue
                        key = row["case"], row["frame"]
                        if key not in local or assisted_rank(row) < assisted_rank(
                            local[key]
                        ):
                            local[key] = row
                summaries[str(path)] = digest(path)
            frames.update(local)
    if not summaries:
        raise ValueError("No complete screen records found")
    provenance = {
        "path_counts": dict(counts),
        "screen_sha256": summaries,
        "source_input_sha256": source_hashes,
    }
    return sorted(frames.values(), key=assisted_rank), provenance, totals


def audit_shortlist(rows: list[dict]) -> list[dict]:
    """Bound native audit to twenty frames while exposing both nearpass classes."""
    chosen = []
    for category, quota in zip(CLASSES, (20, 10, 10, 20), strict=True):
        candidates = [
            row
            for row in rows
            if row["source_class"] == category
            and all(
                row["case"] != prior["case"]
                or abs(row["time_ps"] - prior["time_ps"]) >= 1000
                for prior in chosen
            )
        ]
        chosen.extend(diverse(candidates, min(quota, 20 - len(chosen))))
    return chosen


def select_assisted(rows: list[dict]) -> list[dict]:
    """Select at most two joint candidates and two explicitly diagnostic exports."""
    qualified = sorted(
        [r for r in rows if r["source_class"] == CLASSES[0]],
        key=assisted_rank,
    )
    selected = [
        {**r, "selection_role": "joint_qualified_assisted"}
        for r in diverse(qualified, 2)
    ]
    diagnostics = [
        row
        for row in sorted(rows, key=assisted_rank)
        if row["source_class"] != CLASSES[0]
        and all(
            row["case"] != prior["case"]
            or abs(row["time_ps"] - prior["time_ps"]) >= 1000
            for prior in selected
        )
    ]
    selected.extend(
        {**r, "selection_role": "diagnostic_nearpass_not_qualified"}
        for r in diverse(diagnostics, 2)
    )
    return selected


def review_status(selected: list[dict]) -> str:
    """Report assisted qualification without promoting diagnostic rows."""
    if any(row["selection_role"] == "joint_qualified_assisted" for row in selected):
        return "joint_qualified_assisted_found"
    return "no_joint_qualified_assisted_candidate"


def review_source_hashes(script: Path) -> dict[str, str]:
    """Record the review CLI and relay helper versions used for interpretation."""
    modules = script.resolve().parents[1] / "src/nyl12_covalent"
    return {
        str(path.resolve()): digest(path)
        for path in [script, *modules.glob("leaving_group_relay*.py")]
    }
