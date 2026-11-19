"""Check two-leg hydrogen-bond geometry for candidate Nyl12 relay waters."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_angles, calc_bonds

LOGGER = logging.getLogger(__name__)
PAIRS = ("OG1_ThrN", "OG1_leavingN")
DA_MAX_A = 3.5
HA_MAX_A = 2.5
ANGLE_MIN_DEG = 120.0


def atom_index(universe: mda.Universe, resid: int, name: str) -> int:
    """Map a residue/atom name to exactly one topology index.

    Args:
        universe: Loaded topology.
        resid: One-based Amber residue number.
        name: Amber atom name.

    Returns:
        Zero-based atom index.
    """
    selected = universe.select_atoms(f"resid {resid} and name {name}")
    if len(selected) != 1:
        raise ValueError(f"Expected one {resid}:{name}, got {len(selected)}")
    return int(selected.indices[0])


def hbond(
    xyz: np.ndarray, donor: int, hydrogen: int, acceptor: int, box: np.ndarray
) -> tuple[bool, float, float, float]:
    """Evaluate donor–hydrogen–acceptor geometry with periodic boundaries.

    Args:
        xyz: Current frame coordinates.
        donor: Donor heavy-atom index.
        hydrogen: Attached donor hydrogen index.
        acceptor: Acceptor heavy-atom index.
        box: Current periodic cell.

    Returns:
        Pass flag, donor–acceptor, H–acceptor, and D–H–A angle.
    """
    d, h, a = (xyz[index].reshape(1, 3) for index in (donor, hydrogen, acceptor))
    da = float(calc_bonds(d, a, box=box)[0])
    ha = float(calc_bonds(h, a, box=box)[0])
    angle = float(np.degrees(calc_angles(d, h, a, box=box)[0]))
    return da <= DA_MAX_A and ha <= HA_MAX_A and angle >= ANGLE_MIN_DEG, da, ha, angle


def events_by_frame(path: Path) -> dict[int, list[tuple[str, int, float]]]:
    """Index candidate oxygen contacts by trajectory frame.

    Args:
        path: Oxygen candidate event CSV.

    Returns:
        Frame-indexed pair, water residue and source time records.
    """
    records: dict[int, list[tuple[str, int, float]]] = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["pair"] not in PAIRS:
                raise ValueError(f"Unexpected pair: {row['pair']}")
            records[int(row["frame"])].append(
                (row["pair"], int(row["water_resid"]), float(row["time_ps"]))
            )
    return records


def analyze_case(
    case: dict[str, Any], events: Path, geometry: Path, output: Path
) -> dict[str, Any]:
    """Screen hydrogen-bond direction in all oxygen-candidate frames.

    Args:
        case: Matched MD topology and trajectory manifest entry.
        events: Oxygen candidate event CSV.
        geometry: Per-frame Tier geometry CSV.
        output: New detailed hydrogen-bond CSV.

    Returns:
        Per-pair counts of complete two-leg candidates and Tier A overlap.
    """
    universe = mda.Universe(case["topology"], case["trajectories"])
    if len(universe.trajectory) != 10000:
        raise ValueError("Expected 10,000 production frames")
    og1 = atom_index(universe, 929, "OG1")
    hg1 = atom_index(universe, 929, "HG1")
    targets = {
        "OG1_ThrN": atom_index(universe, 929, "N"),
        "OG1_leavingN": atom_index(universe, 1315, "N"),
    }
    with geometry.open(newline="") as handle:
        tier_a = {
            int(row["frame"]) for row in csv.DictReader(handle) if row["tier"] == "A"
        }
    frame_events = events_by_frame(events)
    water_resids = {
        resid for entries in frame_events.values() for _, resid, _ in entries
    }
    water_atoms: dict[int, tuple[int, int, int]] = {}
    for resid in water_resids:
        residue = universe.residues[resid - 1]
        if residue.resid != resid or residue.resname != "WAT":
            raise ValueError(f"Water residue mapping failed at {resid}")
        atoms = {atom.name: int(atom.index) for atom in residue.atoms}
        if set(atoms) != {"O", "H1", "H2"}:
            raise ValueError(f"Unexpected water inventory at {resid}")
        water_atoms[resid] = atoms["O"], atoms["H1"], atoms["H2"]
    complete: dict[str, set[int]] = {pair: set() for pair in PAIRS}
    complete_ids: dict[str, dict[int, list[int]]] = {
        pair: defaultdict(list) for pair in PAIRS
    }
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            (
                "frame",
                "time_ps",
                "pair",
                "water_resid",
                "og1_to_water_hbond",
                "water_to_n_hbond",
                "complete_relay",
                "og1_h_to_water_A",
                "og1_h_water_angle_deg",
                "best_water_h_to_n_A",
                "best_water_h_n_angle_deg",
            )
        )
        for frame in sorted(frame_events):
            ts = universe.trajectory[frame - 1]
            xyz = ts.positions
            for pair, water_resid, source_time in frame_events[frame]:
                if abs(float(ts.time) - source_time) > 1e-3:
                    raise ValueError("Event and trajectory time disagree")
                water_o, *water_h = water_atoms[water_resid]
                first = hbond(xyz, og1, hg1, water_o, ts.dimensions)
                second = [
                    hbond(xyz, water_o, h, targets[pair], ts.dimensions)
                    for h in water_h
                ]
                best = min(second, key=lambda item: item[2])
                water_to_n = any(item[0] for item in second)
                is_complete = first[0] and water_to_n
                if is_complete:
                    complete[pair].add(frame)
                    complete_ids[pair][water_resid].append(frame)
                writer.writerow(
                    (
                        frame,
                        float(ts.time),
                        pair,
                        water_resid,
                        first[0],
                        water_to_n,
                        is_complete,
                        first[2],
                        first[3],
                        best[2],
                        best[3],
                    )
                )
    longest: dict[str, int] = {}
    for pair in PAIRS:
        longest[pair] = 0
        for frames in complete_ids[pair].values():
            streak = 0
            previous = -2
            for frame in frames:
                streak = streak + 1 if frame == previous + 1 else 1
                longest[pair] = max(longest[pair], streak)
                previous = frame
    return {
        "case": case["name"],
        "tier_a_frames": len(tier_a),
        "hbond_criteria": {
            "donor_acceptor_max_A": DA_MAX_A,
            "hydrogen_acceptor_max_A": HA_MAX_A,
            "donor_hydrogen_acceptor_min_deg": ANGLE_MIN_DEG,
        },
        "complete_two_leg_frames": {pair: len(complete[pair]) for pair in PAIRS},
        "complete_two_leg_tier_a_frames": {
            pair: len(complete[pair] & tier_a) for pair in PAIRS
        },
        "longest_same_water_streak_frames": longest,
        "note": "Geometric hydrogen bonds in fixed-charge MD; no proton-transfer event.",
    }


def main() -> None:
    """Run the three-case hydrogen-bond direction audit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--md-manifest", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--water", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    cases = {x["name"]: x for x in json.loads(args.md_manifest.read_text())["cases"]}
    geometry = {
        x["name"]: x for x in json.loads(args.comparison.read_text())["replicas"]
    }
    water = json.loads((args.water / "summary.json").read_text())
    results = []
    for item in water:
        name = item["case"]
        result = analyze_case(
            cases[name],
            args.water / f"{name}_events.csv",
            Path(geometry[name]["frame_csv"]),
            args.output / f"{name}_hbonds.csv",
        )
        results.append(result)
        LOGGER.info("%s: %s", name, result["complete_two_leg_tier_a_frames"])
    (args.output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
