"""Relate candidate water oxygen contacts to matched Tier A frames."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

PAIRS = ("OG1_ThrN", "OG1_leavingN")


def summarize_case(geometry_csv: Path, events_csv: Path) -> dict[str, object]:
    """Count water candidate frames within and outside catalytic Tier A.

    Args:
        geometry_csv: Matched all-frame Tier geometry audit.
        events_csv: Candidate water oxygen events.

    Returns:
        Frame counts and conditional fractions.
    """
    with geometry_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    frames = {int(row["frame"]) for row in rows}
    tier_a = {int(row["frame"]) for row in rows if row["tier"] == "A"}
    if len(frames) != 10000 or frames != set(range(1, 10001)):
        raise ValueError("Expected 10,000 contiguous production frames")
    events = {pair: set() for pair in PAIRS}
    with events_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            pair = row["pair"]
            if pair not in events:
                raise ValueError(f"Undeclared water pair {pair}")
            frame = int(row["frame"])
            if frame not in frames:
                raise ValueError(f"Water event has invalid frame {frame}")
            events[pair].add(frame)
    return {
        "tier_a_frames": len(tier_a),
        "pairs": {
            pair: {
                "all_candidate_frames": len(contact),
                "tier_a_candidate_frames": len(contact & tier_a),
                "fraction_of_tier_a_with_candidate": (
                    len(contact & tier_a) / len(tier_a) if tier_a else None
                ),
            }
            for pair, contact in events.items()
        },
    }


def main() -> None:
    """Write a per-case conditional water occupancy summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--water", required=True, type=Path)
    args = parser.parse_args()
    summary_path = args.water / "tier_a_water_summary.json"
    if summary_path.exists():
        raise FileExistsError(summary_path)
    comparison = json.loads(args.comparison.read_text())
    by_name = {item["name"]: item for item in comparison["replicas"]}
    water = json.loads((args.water / "summary.json").read_text())
    output = {}
    for case in water:
        name = case["case"]
        output[name] = summarize_case(
            Path(by_name[name]["frame_csv"]), args.water / f"{name}_events.csv"
        )
    summary_path.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
