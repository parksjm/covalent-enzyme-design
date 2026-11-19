#!/usr/bin/env python3
"""Screen modular-PA6 MD frames for native-like catalytic Thr torsions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.modular_es_snapshot_selection import (  # noqa: E402
    classify_geometry,
    read_geometry_timeseries,
)

LOGGER = logging.getLogger(__name__)
LATE_WINDOWS_NS = ((50.0, 66.67), (66.67, 83.34), (83.34, 100.0))


def wrapped_delta(value: float, reference: float) -> float:
    """Return the signed minimum periodic angular difference in degrees."""
    return (value - reference + 180.0) % 360.0 - 180.0


def read_cpptraj_scalar(path: Path) -> list[float]:
    """Read one scalar value per sequential cpptraj output frame."""
    values: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        frame = int(fields[0])
        if frame != len(values) + 1:
            raise ValueError(f"Nonsequential cpptraj frames in {path}")
        values.append(float(fields[1]))
    if not values:
        raise ValueError(f"No scalar values in {path}")
    return values


def cpptraj_text(topology: Path, trajectory: Path, frames: int, output: Path) -> str:
    """Return cpptraj input for Thr233 psi and chi1 time series."""
    if frames < 1:
        raise ValueError("frames must be positive")
    return f"""parm {topology.resolve()}
trajin {trajectory.resolve()} 1 {frames} 1
dihedral psi :929@N :929@CA :929@C :930@N out {output / "psi.dat"}
dihedral chi1 :929@N :929@CA :929@CB :929@OG1 out {output / "chi1.dat"}
run
"""


def representative_joint_frames(
    rows: list[dict[str, Any]],
    windows: tuple[tuple[float, float], ...] = LATE_WINDOWS_NS,
) -> list[dict[str, Any]]:
    """Choose midpoint-representative native-torsion Tier A/B frames."""
    selected = []
    for start, end in windows:
        candidates = [
            row
            for row in rows
            if start <= row["time_ns"] < end
            and row["native_torsions"]
            and row["tier"] in {"A", "B"}
        ]
        if not candidates:
            continue
        tier = "A" if any(row["tier"] == "A" for row in candidates) else "B"
        eligible = [row for row in candidates if row["tier"] == tier]
        midpoint = (start + end) / 2.0
        chosen = min(
            eligible,
            key=lambda row: (abs(row["time_ns"] - midpoint), row["frame"]),
        )
        selected.append(
            {
                "frame": int(chosen["frame"]),
                "time_ns": chosen["time_ns"],
                "tier": tier,
                "window_start_ns": start,
                "window_end_ns": end,
                "joint_eligible_count": len(eligible),
                "psi_delta": chosen["psi_delta"],
                "chi1_delta": chosen["chi1_delta"],
            }
        )
    return selected


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def screen(
    topology: Path,
    trajectory: Path,
    geometry_csv: Path,
    cpptraj: Path,
    output: Path,
    reference_psi: float,
    reference_chi1: float,
) -> dict[str, Any]:
    """Calculate torsions, join geometry, and report native-like Tier A/B frames."""
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    geometry = read_geometry_timeseries(geometry_csv)
    cpptraj_input = output / "torsions.cpptraj"
    cpptraj_input.write_text(
        cpptraj_text(topology, trajectory, len(geometry), output), encoding="utf-8"
    )
    environment = dict(os.environ)
    amberhome = cpptraj.parent.parent
    environment["AMBERHOME"] = str(amberhome)
    environment["OMP_NUM_THREADS"] = "1"
    library = str(amberhome / "lib")
    environment["LD_LIBRARY_PATH"] = (
        f"{library}:{environment['LD_LIBRARY_PATH']}"
        if environment.get("LD_LIBRARY_PATH")
        else library
    )
    with (output / "cpptraj.log").open("w", encoding="utf-8") as log:
        subprocess.run(
            [str(cpptraj), "-i", str(cpptraj_input)],
            check=True,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    psi = read_cpptraj_scalar(output / "psi.dat")
    chi1 = read_cpptraj_scalar(output / "chi1.dat")
    if len(psi) != len(geometry) or len(chi1) != len(geometry):
        raise ValueError("Torsion and geometry frame counts differ")
    rows = []
    for metrics, psi_value, chi1_value in zip(geometry, psi, chi1, strict=True):
        psi_delta = wrapped_delta(psi_value, reference_psi)
        chi1_delta = wrapped_delta(chi1_value, reference_chi1)
        tier = classify_geometry(metrics)
        rows.append(
            {
                **metrics,
                "tier": tier or "none",
                "thr233_psi": psi_value,
                "thr233_chi1": chi1_value,
                "psi_delta": psi_delta,
                "chi1_delta": chi1_delta,
                "native_torsions": abs(psi_delta) <= 35.0 and abs(chi1_delta) <= 20.0,
            }
        )
    with (output / "geometry_and_thr_torsions.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    eligible = [row for row in rows if row["tier"] in {"A", "B"}]
    joint = [row for row in eligible if row["native_torsions"]]
    summary = {
        "schema_version": 1,
        "frame_count": len(rows),
        "reference_degrees": {"psi": reference_psi, "chi1": reference_chi1},
        "native_limits_degrees": {
            "absolute_psi_delta": 35.0,
            "absolute_chi1_delta": 20.0,
        },
        "tier_counts": {
            tier: sum(row["tier"] == tier for row in rows) for tier in ("A", "B")
        },
        "native_torsion_count": sum(row["native_torsions"] for row in rows),
        "tier_a_or_b_and_native_torsion_count": len(joint),
        "joint_frames": [int(row["frame"]) for row in joint],
        "representative_late_joint_frames": representative_joint_frames(rows),
        "sources": {
            str(path.resolve()): _sha256(path)
            for path in (topology, trajectory, geometry_csv, cpptraj, Path(__file__))
        },
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    """Parse CLI arguments and run the torsion screen."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--geometry-csv", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-psi", type=float, default=134.11027642801542)
    parser.add_argument("--reference-chi1", type=float, default=-58.96441787770226)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    report = screen(
        args.topology,
        args.trajectory,
        args.geometry_csv,
        args.cpptraj,
        args.output,
        args.reference_psi,
        args.reference_chi1,
    )
    LOGGER.info(
        "Joint Tier A/B and native-torsion frames: %d",
        report["tier_a_or_b_and_native_torsion_count"],
    )


if __name__ == "__main__":
    main()
