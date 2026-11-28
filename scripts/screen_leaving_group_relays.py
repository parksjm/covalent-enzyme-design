"""Stream existing MD into bounded explicit-hydrogen leaving-amide relay records."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import logging
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
import parmed
from scipy.io import netcdf_file

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nyl12_covalent.leaving_group_relay_mapping import resolve_roles  # noqa: E402
from nyl12_covalent.leaving_group_relays import (  # noqa: E402
    catalytic_geometry,
    screen_paths,
)
from nyl12_covalent.nalpha_water_relay import cell_matrix  # noqa: E402

LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return SHA-256 of a source file."""
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def scan_case(
    case: dict, criteria: dict, output: Path, stride: int, frames: set[int] | None
) -> dict:
    """Stream a case, preserving strict and nearest failed-route evidence.

    Args:
        case: Explicit independent topology, role and geometry mappings.
        criteria: D-A/H-A/angle thresholds.
        output: New output directory.
        stride: Saved-frame sampling stride.
        frames: Optional one-based global saved-frame refinement selection.

    Returns:
        Case summary including immutable input hashes and selected frame list.
    """
    topology = parmed.load_file(case["topology"])
    roles, geometry = resolve_roles(topology, case)
    counts: Counter = Counter()
    snapshots = []
    inputs = [Path(case["topology"]), *map(Path, case["trajectories"])]
    source_stats = {
        str(path): (path.stat().st_size, path.stat().st_mtime_ns) for path in inputs
    }
    hash_started = time.monotonic()
    hashes = {str(path): digest(path) for path in inputs}
    hash_seconds = time.monotonic() - hash_started
    LOGGER.info(
        "%s input hashing completed in %.1f seconds", case["name"], hash_seconds
    )
    csv_path = output / f"{case['name']}_paths.csv"
    json_path = output / f"{case['name']}_paths.jsonl"
    columns = [
        "case",
        "frame",
        "local_frame",
        "trajectory",
        "time_ps",
        "start",
        "strict",
        "water_count",
        "internal_count",
        "violation",
        "d1",
        "d2",
        "d3",
        "attack",
        "projected_phi",
        "path",
        "states",
        "conditional_roles",
        "atom_indices1",
        "edges",
        "leaving_N_acceptor_conditional",
    ]
    with csv_path.open("x", newline="") as ch, json_path.open("x") as jh:
        writer = csv.DictWriter(ch, fieldnames=columns)
        writer.writeheader()
        offset = 0
        for trajectory in case["trajectories"]:
            with netcdf_file(trajectory, mmap=True) as source:
                if source.variables["coordinates"].shape[1:] != (
                    len(topology.atoms),
                    3,
                ):
                    raise ValueError("Trajectory atom count differs from topology")
                count = source.variables["coordinates"].shape[0]
                for local in range(count):
                    frame = offset + local + 1
                    if frames is not None and frame not in frames:
                        continue
                    if frames is None and (frame - 1) % stride:
                        continue
                    xyz = np.array(source.variables["coordinates"][local], dtype=float)
                    cell = cell_matrix(
                        np.array(source.variables["cell_lengths"][local]),
                        np.array(source.variables["cell_angles"][local]),
                    )
                    metrics = catalytic_geometry(xyz, geometry, cell)
                    paths = screen_paths(xyz, roles, geometry["nleave"], cell, criteria)
                    counts["sampled_frames"] += 1
                    strict_frame = False
                    for route in paths:
                        row = {
                            "case": case["name"],
                            "frame": frame,
                            "local_frame": local + 1,
                            "trajectory": trajectory,
                            "time_ps": float(source.variables["time"][local]),
                            **route,
                            **metrics,
                        }
                        jh.write(json.dumps(row, allow_nan=False) + "\n")
                        writer.writerow(
                            {
                                k: json.dumps(v) if isinstance(v, list) else v
                                for k, v in row.items()
                            }
                        )
                        counts[
                            "strict_paths" if route["strict"] else "nearest_nearpaths"
                        ] += 1
                        strict_frame |= route["strict"]
                    if strict_frame:
                        counts["strict_frames"] += 1
                        snapshots.append(
                            {
                                "case": case["name"],
                                "frame": frame,
                                "trajectory": trajectory,
                                "local_frame": local + 1,
                            }
                        )
                    if counts["sampled_frames"] % 100 == 0:
                        LOGGER.info("%s frame %d %s", case["name"], frame, dict(counts))
                        ch.flush()
                        jh.flush()
                offset += count
    if frames is not None and any(frame < 1 or frame > offset for frame in frames):
        raise ValueError("Requested frame outside trajectory inventory")
    if any(
        (path.stat().st_size, path.stat().st_mtime_ns) != source_stats[str(path)]
        for path in inputs
    ):
        raise ValueError("Source input changed during screen")
    water_inventory = [r for r in roles if r["water"]]
    return {
        "hash_seconds": hash_seconds,
        "input_size_mtime_ns": source_stats,
        "water_count": len(water_inventory),
        "water_inventory_sha256": hashlib.sha256(
            json.dumps(water_inventory, sort_keys=True).encode()
        ).hexdigest(),
        "case": case["name"],
        "counts": dict(counts),
        "total_frames": offset,
        "input_sha256": hashes,
        "selected_snapshots": snapshots,
        "role_mapping": [role for role in roles if not role["water"]],
        "geometry_mapping": geometry,
    }


def main() -> None:
    """Run a new immutable screen from an explicit JSON configuration."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stride", type=int)
    parser.add_argument(
        "--frames", type=Path, help="JSON list of {case,frame} for refinement"
    )
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    config = json.loads(args.config.read_text())
    stride = args.stride if args.stride is not None else config["stride"]
    if stride < 1:
        raise ValueError("Stride must be positive")
    criteria = config["criteria"]
    if (
        criteria["da_A"] <= 0
        or criteria["ha_A"] <= 0
        or not 0 < criteria["angle_deg"] <= 180
    ):
        raise ValueError("Invalid explicit hydrogen-bond thresholds")
    selected = json.loads(args.frames.read_text()) if args.frames else None
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    summary: dict[str, Any] = {
        "criteria": criteria,
        "stride": stride,
        "cases": [],
        "interpretation": "Structural reaction-preparation hypotheses; "
        "conditional leaving-amide-N acceptor; no validated barrier.",
    }
    for case in config["cases"]:
        frames = (
            {int(row["frame"]) for row in selected if row["case"] == case["name"]}
            if selected is not None
            else None
        )
        result = scan_case(case, criteria, args.output, stride, frames)
        summary["cases"].append(result)
        (args.output / "progress.json").write_text(json.dumps(summary, indent=2) + "\n")
        LOGGER.info("Completed %s: %s", case["name"], result["counts"])
    sources = [
        Path(__file__),
        *Path(__file__)
        .resolve()
        .parents[1]
        .joinpath("src/nyl12_covalent")
        .glob("leaving_group_relay*.py"),
    ]
    summary["source_sha256"] = {str(p): digest(p) for p in sources}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    snapshots = [row for case in summary["cases"] for row in case["selected_snapshots"]]
    (args.output / "selected_snapshots.json").write_text(
        json.dumps(snapshots, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
