"""Measure primary nearpasses and independent native-gate failure distributions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import sys

import numpy as np
import parmed
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from nyl12_covalent.nalpha_water_relay import (  # noqa:E402
    angle,
    cell_matrix,
    hbond,
    minimum_image,
)
from screen_nalpha_leaving_water import digest, mapping  # noqa:E402

LOGGER = logging.getLogger(__name__)


def summarize_failures(directory: Path) -> dict:
    """Count overlapping native-gate failures in already screened Tier A/B frames."""
    summary = {}
    for path in sorted(directory.glob("*_gates.csv")):
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        checks = {
            "OG1_shift_gt075A": ("native_og1_shift_A", 0.75),
            "Nalpha_shift_gt125A": ("native_nalpha_shift_A", 1.25),
            "abs_psi_delta_gt35deg": ("psi_delta_deg", 35),
            "abs_chi1_delta_gt20deg": ("chi1_delta_deg", 20),
        }
        summary[path.stem.removesuffix("_gates")] = {
            "tier_ab_positive_phi_frames": len(rows),
            "native_pass_frames": sum(r["native_pass"] == "True" for r in rows),
            "failure_counts_overlap": {
                label: sum(abs(float(r[column])) > limit for r in rows)
                for label, (column, limit) in checks.items()
            },
            "median_metrics": {
                column: float(np.median([float(r[column]) for r in rows]))
                for column, _ in checks.values()
            }
            if rows
            else {},
            "source": str(path),
            "sha256": digest(path),
        }
    return summary


def primary_nearpasses(config: dict, primary: Path) -> list[dict]:
    """Report the nearest observed water bridge for every fully eligible primary frame."""
    records = []
    for case in config["cases"]:
        topology = parmed.load_file(case["topology"])
        mapped = mapping(topology, config)
        waters = np.asarray(mapped["waters"], dtype=int)
        with (primary / f"{case['name']}_gates.csv").open(newline="") as handle:
            frames = {
                int(r["frame"]): r
                for r in csv.DictReader(handle)
                if r["native_pass"] == "True" and r["intact_pa6"] == "True"
            }
        with Path(case["geometry_csv"]).open(newline="") as handle:
            geometry = {int(r["frame"]): r for r in csv.DictReader(handle)}
        offset = 0
        for trajectory in case["trajectories"]:
            with netcdf_file(trajectory, mmap=True) as source:
                count = source.variables["coordinates"].shape[0]
                for frame in sorted(frames):
                    local = frame - 1 - offset
                    if not 0 <= local < count:
                        continue
                    xyz = np.array(source.variables["coordinates"][local], dtype=float)
                    cell = cell_matrix(
                        np.array(source.variables["cell_lengths"][local]),
                        np.array(source.variables["cell_angles"][local]),
                    )
                    da = np.linalg.norm(
                        minimum_image(xyz[waters[:, 1]] - xyz[mapped["nalpha"]], cell),
                        axis=1,
                    )
                    dl = np.linalg.norm(
                        minimum_image(xyz[waters[:, 1]] - xyz[mapped["nleave"]], cell),
                        axis=1,
                    )
                    best = int(np.argmin(np.maximum(da, dl)))
                    resid, oxygen, h1, h2 = waters[best]
                    legs = [
                        hbond(xyz, oxygen, h, mapped["nleave"], cell) for h in (h1, h2)
                    ]
                    outgoing = min(legs, key=lambda r: (not r["pass"], r["ha_A"]))
                    incoming = min(
                        [
                            hbond(xyz, mapped["nalpha"], h, oxygen, cell)
                            for h in mapped["nalpha_h"]
                        ],
                        key=lambda r: (not r["pass"], r["ha_A"]),
                    )
                    local_indices = [
                        value for value in mapped.values() if isinstance(value, int)
                    ]
                    local_points = minimum_image(
                        xyz[local_indices] - xyz[mapped["carbon"]], cell
                    )
                    points = dict(zip(local_indices, local_points, strict=True))
                    records.append(
                        {
                            "case": case["name"],
                            "frame": frame,
                            "time_ps": float(geometry[frame]["time_ps"]),
                            "water_resid": int(resid),
                            "nalpha_nleave_A": float(
                                np.linalg.norm(
                                    minimum_image(
                                        xyz[mapped["nalpha"]] - xyz[mapped["nleave"]],
                                        cell,
                                    )
                                )
                            ),
                            "nalpha_water_A": float(da[best]),
                            "water_nleave_A": float(dl[best]),
                            "max_bridge_leg_A": float(max(da[best], dl[best])),
                            "water_h_to_nleave_A": outgoing["ha_A"],
                            "water_h_to_nleave_angle_deg": outgoing["angle_deg"],
                            "nalpha_h_to_water_A": incoming["ha_A"],
                            "nalpha_h_to_water_angle_deg": incoming["angle_deg"],
                            "water_to_nleave_hbond": outgoing["pass"],
                            "nalpha_to_water_hbond": incoming["pass"],
                            **{
                                k: float(geometry[frame][k])
                                for k in ("d1", "d2", "d3", "attack", "projected_phi")
                            },
                            "cb_og1_carbonyl_angle_deg": angle(
                                *[points[mapped[k]] for k in ("cb", "og1", "carbon")]
                            ),
                            "nleave_scissile_c_c1_angle_deg": angle(
                                *[
                                    points[mapped[k]]
                                    for k in ("carbon", "nleave", "leaving_c1")
                                ]
                            ),
                            "native_thr_rmsd_A": float(
                                frames[frame]["native_thr_rmsd_A"]
                            ),
                            "diagnostic_only": True,
                        }
                    )
                offset += count
    return sorted(records, key=lambda r: (r["max_bridge_leg_A"], r["frame"]))


def main() -> None:
    """Write diagnostics separately without extracting or promoting failed frames."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--controls", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    config = json.loads(args.config.read_text())
    records = primary_nearpasses(config, args.primary)
    with (args.output / "primary_nearpasses.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(records[0]) if records else ["case", "frame"]
        )
        writer.writeheader()
        writer.writerows(records)
    result = {
        "primary": summarize_failures(args.primary),
        "controls": summarize_failures(args.controls),
        "primary_nearpasses": records,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "acceptance_unchanged": True,
        "diagnostic_only": True,
    }
    (args.output / "diagnostics.json").write_text(json.dumps(result, indent=2) + "\n")
    LOGGER.info("Primary nearpass diagnostics: %d frames", len(records))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
