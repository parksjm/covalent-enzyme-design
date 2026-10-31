#!/usr/bin/env python3
"""Analyze productive geometry in a fixed prefix of a modular PA6 trajectory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import statistics
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.active_site_geometry import projected_herschlag_phi  # noqa: E402

LOGGER = logging.getLogger("analyze_modular_pa6_md_trajectory")


@dataclass(frozen=True)
class Thresholds:
    """Five-metric operational geometry thresholds."""

    d1: float
    oxy: float
    angle_min: float
    angle_max: float
    phi_min: float
    phi_max: float


TIERS = {
    "tier_a": Thresholds(3.2, 3.8, 90.0, 120.0, 70.0, 110.0),
    "tier_b": Thresholds(3.5, 4.0, 85.0, 125.0, 60.0, 125.0),
    "tier_c": Thresholds(3.8, 4.2, 80.0, 130.0, 45.0, 140.0),
}


def cpptraj_text(topology: Path, trajectory: Path, output: Path, frames: int) -> str:
    """Return cpptraj input for a stable trajectory prefix and five metrics."""
    selection = "(:929@OG1,:1314@C,:1314@O,:1315@N,:1314@C5)"
    return f"""parm {topology.resolve()}
trajin {trajectory.resolve()} 1 {frames} 1
autoimage anchor :1313-1316
rms backbone :1-1312@CA,C,N,O first out {output / "backbone_rmsd.dat"}
rms ligand :1313-1316&!@H= first nofit out {output / "ligand_positional_rmsd.dat"}
distance d1 :929@OG1 :1314@C out {output / "d1.dat"}
distance d2 :113@N :1314@O out {output / "d2.dat"}
distance d3 :186@ND2 :1314@O out {output / "d3.dat"}
angle attack :929@OG1 :1314@C :1314@O out {output / "attack.dat"}
trajout {output / "prefix_35ns.nc"} netcdf
strip !{selection}
trajout {output / "projected_phi_atoms.pdb"} pdb multi nobox
run
"""


def read_dat(path: Path) -> dict[int, float]:
    """Read a cpptraj scalar dataset keyed by one-based output frame."""
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if fields and fields[0].isdigit():
            values[int(fields[0])] = float(fields[1])
    if not values:
        raise ValueError(f"No data rows in {path}")
    return values


def read_projected_phi(path: Path) -> dict[int, float]:
    """Calculate projected Herschlag phi from cpptraj five-atom PDB frames."""
    values = {}
    split_paths = sorted(
        path.parent.glob(f"{path.name}.*"),
        key=lambda item: int(item.name.rsplit(".", 1)[-1]),
    )
    requirements = [
        (929, "OG1"),
        (1314, "O"),
        (1314, "C"),
        (1315, "N"),
        (1314, "C5"),
    ]
    for frame_path in split_paths:
        coordinates = {}
        for line in frame_path.read_text(encoding="ascii").splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                key = (int(line[22:26]), line[12:16].strip())
                coordinates[key] = tuple(
                    float(line[start:end])
                    for start, end in ((30, 38), (38, 46), (46, 54))
                )
        if set(coordinates) != set(requirements):
            raise ValueError(f"Projected-phi atom mismatch in {frame_path}")
        points = [coordinates[key] for key in requirements]
        frame = int(frame_path.name.rsplit(".", 1)[-1])
        values[frame] = projected_herschlag_phi(*points)
    if not values:
        raise ValueError("No projected-phi coordinate frames were written")
    return values


def passes(row: dict[str, float], thresholds: Thresholds) -> bool:
    """Return whether one frame passes all five geometry thresholds."""
    return (
        row["d1"] <= thresholds.d1
        and row["d2"] <= thresholds.oxy
        and row["d3"] <= thresholds.oxy
        and thresholds.angle_min <= row["attack"] <= thresholds.angle_max
        and thresholds.phi_min <= abs(row["projected_phi"]) <= thresholds.phi_max
    )


def summarize(rows: list[dict[str, float]]) -> dict[str, Any]:
    """Summarize geometry occupancy, drift, and positional ligand RMSD."""
    midpoint = len(rows) // 2
    summary: dict[str, Any] = {
        "n_frames": len(rows),
        "last_time_ns": rows[-1]["time_ns"],
        "mean_d1_A": statistics.fmean(row["d1"] for row in rows),
        "mean_d2_A": statistics.fmean(row["d2"] for row in rows),
        "mean_d3_A": statistics.fmean(row["d3"] for row in rows),
        "mean_attack_deg": statistics.fmean(row["attack"] for row in rows),
        "median_abs_projected_phi_deg": statistics.median(
            abs(row["projected_phi"]) for row in rows
        ),
        "median_ligand_positional_rmsd_A": statistics.median(
            row["ligand_rmsd"] for row in rows
        ),
        "final_ligand_positional_rmsd_A": rows[-1]["ligand_rmsd"],
    }
    for name, thresholds in TIERS.items():
        summary[f"{name}_fraction"] = statistics.fmean(
            passes(row, thresholds) for row in rows
        )
        summary[f"{name}_first_half"] = statistics.fmean(
            passes(row, thresholds) for row in rows[:midpoint]
        )
        summary[f"{name}_second_half"] = statistics.fmean(
            passes(row, thresholds) for row in rows[midpoint:]
        )
    summary["interim_persistence_rule_passed"] = (
        summary["tier_b_fraction"] >= 0.20
        and summary["tier_c_fraction"] >= 0.50
        and summary["mean_d1_A"] <= 4.0
    )
    return summary


def analyze(
    topology: Path,
    trajectory: Path,
    output: Path,
    cpptraj: Path,
    variant: str,
    replica: str,
    frames: int,
) -> dict[str, Any]:
    """Run cpptraj and write an immutable per-frame and summary report."""
    if output.exists():
        raise FileExistsError(output)
    if frames < 2:
        raise ValueError("At least two frames are required")
    output.mkdir(parents=True)
    source_stat = trajectory.stat()
    script = output / "analysis.cpptraj"
    script.write_text(cpptraj_text(topology, trajectory, output, frames))
    environment = dict(os.environ)
    environment["AMBERHOME"] = str(cpptraj.parents[1])
    amber_lib = str(cpptraj.parents[1] / "lib")
    environment["LD_LIBRARY_PATH"] = (
        f"{amber_lib}:{environment.get('LD_LIBRARY_PATH', '')}"
    )
    with (output / "cpptraj.log").open("w") as log:
        subprocess.run(
            [str(cpptraj), "-i", str(script)],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
        )
    datasets = {
        name: read_dat(output / filename)
        for name, filename in {
            "d1": "d1.dat",
            "d2": "d2.dat",
            "d3": "d3.dat",
            "attack": "attack.dat",
            "ligand_rmsd": "ligand_positional_rmsd.dat",
        }.items()
    }
    datasets["projected_phi"] = read_projected_phi(output / "projected_phi_atoms.pdb")
    frame_ids = sorted(set.intersection(*(set(values) for values in datasets.values())))
    if len(frame_ids) != frames:
        raise ValueError(f"Expected {frames} complete frames, found {len(frame_ids)}")
    rows = [
        {
            "frame": float(frame),
            "time_ns": (frame - 1) * 0.05,
            **{name: values[frame] for name, values in datasets.items()},
        }
        for frame in frame_ids
    ]
    with (output / "geometry_timeseries.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema_version": 1,
        "status": "interim_fixed_prefix_not_final_100ns",
        "variant": variant,
        "replica": replica,
        "source_trajectory": str(trajectory.resolve()),
        "source_size_at_start": source_stat.st_size,
        "requested_prefix_frames": frames,
        "frame_spacing_ns": 0.05,
        "topology_sha256": _sha256(topology),
        "prefix_trajectory_sha256": _sha256(output / "prefix_35ns.nc"),
        "thresholds": {name: asdict(value) for name, value in TIERS.items()},
        "summary": summarize(rows),
    }
    (output / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    return report


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    """Parse arguments and analyze one trajectory prefix."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--replica", required=True)
    parser.add_argument("--frames", type=int, default=700)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    report = analyze(
        args.topology,
        args.trajectory,
        args.output,
        args.cpptraj,
        args.variant,
        args.replica,
        args.frames,
    )
    LOGGER.info("wrote %s frames to %s", report["requested_prefix_frames"], args.output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileExistsError,
        FileNotFoundError,
        subprocess.CalledProcessError,
        ValueError,
    ) as exc:
        LOGGER.error("%s", exc)
        raise SystemExit(2) from exc
