"""Screen candidate relay-water oxygen geometry in matched MD replicas."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array

LOGGER = logging.getLogger(__name__)
CUTOFF_A = 3.5
PAIRS = ("OG1_ThrN", "OG1_leavingN")


def one_atom(universe: mda.Universe, selection: str) -> mda.AtomGroup:
    """Return exactly one mapped active-site atom.

    Args:
        universe: MD topology and trajectory.
        selection: MDAnalysis atom selection.

    Returns:
        One-atom group.
    """
    group = universe.select_atoms(selection)
    if len(group) != 1:
        raise ValueError(f"Expected one atom for {selection}, found {len(group)}")
    return group


def analyze_case(case: dict[str, object], output: Path) -> dict[str, object]:
    """Count candidate water oxygens and uninterrupted identity streaks.

    Args:
        case: Declared topology and two production trajectories.
        output: New event CSV path.

    Returns:
        Frame occupancy and maximum sampled water-identity streak.
    """
    universe = mda.Universe(str(case["topology"]), case["trajectories"])
    expected = int(case["expected_frames"])
    if len(universe.trajectory) != expected:
        raise ValueError("Unexpected trajectory frame count")
    anchors = (
        one_atom(universe, "resid 929 and name OG1"),
        one_atom(universe, "resid 929 and name N"),
        one_atom(universe, "resid 1315 and name N"),
    )
    waters = universe.select_atoms("resname WAT and name O")
    if len(waters) < 100:
        raise ValueError("Water oxygen inventory is unexpectedly small")
    present_frames = dict.fromkeys(PAIRS, 0)
    longest = dict.fromkeys(PAIRS, 0)
    streaks: dict[str, dict[int, int]] = {name: {} for name in PAIRS}
    previous_time: float | None = None
    cadence: float | None = None
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("frame", "time_ps", "pair", "water_resid", "d1_A", "d2_A"))
        for frame, ts in enumerate(universe.trajectory, start=1):
            if previous_time is not None:
                interval = float(ts.time - previous_time)
                if cadence is None:
                    cadence = interval
                elif abs(interval - cadence) > 1e-3:
                    raise ValueError("Nonuniform water-screen frame cadence")
            previous_time = float(ts.time)
            centers = np.vstack([group.positions[0] for group in anchors])
            distances = distance_array(centers, waters.positions, box=ts.dimensions)
            for name, second in zip(PAIRS, (1, 2), strict=True):
                matches = np.flatnonzero(
                    (distances[0] <= CUTOFF_A) & (distances[second] <= CUTOFF_A)
                )
                if len(matches):
                    present_frames[name] += 1
                current: dict[int, int] = {}
                for local in matches:
                    resid = int(waters[int(local)].resid)
                    current[resid] = streaks[name].get(resid, 0) + 1
                    longest[name] = max(longest[name], current[resid])
                    writer.writerow(
                        (
                            frame,
                            float(ts.time),
                            name,
                            resid,
                            float(distances[0, local]),
                            float(distances[second, local]),
                        )
                    )
                streaks[name] = current
    return {
        "case": case["name"],
        "frames": expected,
        "water_oxygen_cutoff_A": CUTOFF_A,
        "candidate_pair_occupancy": {
            name: present_frames[name] / expected for name in PAIRS
        },
        "longest_identity_streak_frames": longest,
        "frame_cadence_ps": cadence,
        "note": "Heavy-atom proximity only; not a hydrogen bond or proton-transfer path.",
    }


def main() -> None:
    """Screen the three declared pKa snapshot source replicas."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--md-manifest", type=Path, required=True)
    parser.add_argument("--screen-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    cases = {
        case["name"]: case for case in json.loads(args.md_manifest.read_text())["cases"]
    }
    names = [run["case"] for run in json.loads(args.screen_manifest.read_text())]
    if len(names) != 3 or len(set(names)) != 3:
        raise ValueError("Expected three unique declared cases")
    summaries = []
    for name in names:
        result = analyze_case(cases[name], args.output / f"{name}_events.csv")
        summaries.append(result)
        LOGGER.info("%s: %s", name, result["candidate_pair_occupancy"])
    (args.output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
