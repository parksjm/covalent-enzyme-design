"""Per-frame water and reactive-site geometry measurements."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.deacylation_water_geometry import (
    ASP_RESIDUES,
    ATTACK_C,
    CARBONYL_O,
    LEAVING_N,
    NALPHA,
    OG1,
    _hbond,
    angle,
    image_near,
)

PER_FRAME_FIELDS = [
    "frame_index",
    "stage",
    "time_ps",
    "water_residue_index_1based",
    "water_oxygen_atom_index_1based",
    "is_frame_nearest",
    "is_near_attack",
    "O_C_A",
    "O_C_Ocarbonyl_angle_deg",
    "O_Nalpha_A",
    "O_OG1_A",
    "O_Nleaving_A",
    "C_OG1_A",
    "C_Nleaving_A",
    "Nalpha_NH_bond_lengths_A",
    "Nleaving_NH_bond_lengths_A",
    "waterO_Nalpha_D_A",
    "waterH_Nalpha_A_A",
    "waterO_H_Nalpha_angle_deg",
    "waterH_Nalpha_closest_acceptor_atom_index_1based",
    "waterO_Asp272_D_A",
    "waterH_Asp272_A_A",
    "waterO_H_Asp272_angle_deg",
    "waterH_Asp272_closest_acceptor_atom_index_1based",
    "waterO_Asp274_D_A",
    "waterH_Asp274_A_A",
    "waterO_H_Asp274_angle_deg",
    "waterH_Asp274_closest_acceptor_atom_index_1based",
    "O_Asp968_nearestO_A",
    "O_Asp970_nearestO_A",
]


def measure_frames(
    frames: list[dict[str, Any]],
    waters: list[dict[str, Any]],
    topology: Any,
    nalpha_h: list[int],
    nleave_h: list[int],
    output: Path,
) -> tuple[list[tuple[Any, ...]], dict[int, dict[str, Any]]]:
    """Measure nearest waters, all-water minima, and descriptive H-bond geometry."""
    per_frame_path = output / "per_frame_metrics.csv"
    minima: dict[int, dict[str, Any]] = {
        water["o"]: {
            "min_attack_distance_A": math.inf,
            "min_attack_angle_deg": None,
            "near_attack_frame_count": 0,
            "first_near_attack_time_ps": None,
            "stages": {},
        }
        for water in waters
    }
    nearest_rows: list[tuple[Any, ...]] = []
    with per_frame_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=PER_FRAME_FIELDS)
        writer.writeheader()
        for frame_idx, frame in enumerate(frames):
            xyz, box = frame["xyz"], frame["box"]
            carbon = xyz[ATTACK_C - 1]
            carbonyl_o = image_near(xyz[[CARBONYL_O - 1]], carbon, box)[0]
            water_xyz = image_near(xyz[[water["o"] for water in waters]], carbon, box)
            distances = np.linalg.norm(water_xyz - carbon, axis=1)
            angles = angle(
                water_xyz, carbon, np.broadcast_to(carbonyl_o, water_xyz.shape)
            )
            nearest_index = int(np.argmin(distances))
            near_mask = (distances <= 3.5) & (angles >= 95) & (angles <= 125)
            for index, water in enumerate(waters):
                item = minima[water["o"]]
                if distances[index] < item.get("_min_distance", math.inf):
                    item["_min_distance"] = float(distances[index])
                    item["min_attack_distance_A"] = float(distances[index])
                    item["min_attack_angle_deg"] = float(angles[index])
                    item["min_distance_frame_index"] = frame_idx
                    item["min_distance_time_ps"] = frame["time_ps"]
                near = bool(near_mask[index])
                item["near_attack_frame_count"] += int(near)
                stage_counts = item["stages"]
                stage_counts[frame["stage"]] = stage_counts.get(
                    frame["stage"], 0
                ) + int(near)
                if near and item["first_near_attack_time_ps"] is None:
                    item["first_near_attack_time_ps"] = frame["time_ps"]
            nearest = waters[nearest_index]
            nearest_rows.append(
                (
                    frame_idx,
                    frame["stage"],
                    frame["time_ps"],
                    nearest["residue"],
                    nearest["o"] + 1,
                    distances[nearest_index],
                    angles[nearest_index],
                )
            )
            selected = [nearest_index] + [
                index
                for index in range(len(waters))
                if index != nearest_index and near_mask[index]
            ]
            for index in selected:
                water = waters[index]
                oxygen = water_xyz[index]
                hydrogens = image_near(xyz[water["h"]], oxygen, box)
                acceptors = {
                    "Nalpha": ([NALPHA - 1], [NALPHA]),
                    "Asp272": (
                        [atom - 1 for atom in ASP_RESIDUES[968]],
                        list(ASP_RESIDUES[968]),
                    ),
                    "Asp274": (
                        [atom - 1 for atom in ASP_RESIDUES[970]],
                        list(ASP_RESIDUES[970]),
                    ),
                }
                hbond = {
                    name: _hbond(oxygen, hydrogens, xyz[indices], atom_ids, box)
                    for name, (indices, atom_ids) in acceptors.items()
                }
                nalpha_nh = {
                    topology.atoms[atom].name: float(
                        np.linalg.norm(
                            image_near(xyz[[atom]], xyz[NALPHA - 1], box)[0]
                            - xyz[NALPHA - 1]
                        )
                    )
                    for atom in nalpha_h
                }
                nleave_nh = {
                    topology.atoms[atom].name: float(
                        np.linalg.norm(
                            image_near(xyz[[atom]], xyz[LEAVING_N - 1], box)[0]
                            - xyz[LEAVING_N - 1]
                        )
                    )
                    for atom in nleave_h
                }
                carbon_og1 = float(
                    np.linalg.norm(image_near(xyz[[OG1 - 1]], carbon, box)[0] - carbon)
                )
                carbon_nleave = float(
                    np.linalg.norm(
                        image_near(xyz[[LEAVING_N - 1]], carbon, box)[0] - carbon
                    )
                )
                row = {
                    "frame_index": frame_idx,
                    "stage": frame["stage"],
                    "time_ps": frame["time_ps"],
                    "water_residue_index_1based": water["residue"],
                    "water_oxygen_atom_index_1based": water["o"] + 1,
                    "is_frame_nearest": index == nearest_index,
                    "is_near_attack": bool(near_mask[index]),
                    "O_C_A": float(distances[index]),
                    "O_C_Ocarbonyl_angle_deg": float(angles[index]),
                    "O_Nalpha_A": float(
                        np.linalg.norm(
                            image_near(xyz[[NALPHA - 1]], oxygen, box)[0] - oxygen
                        )
                    ),
                    "O_OG1_A": float(
                        np.linalg.norm(
                            image_near(xyz[[OG1 - 1]], oxygen, box)[0] - oxygen
                        )
                    ),
                    "O_Nleaving_A": float(
                        np.linalg.norm(
                            image_near(xyz[[LEAVING_N - 1]], oxygen, box)[0] - oxygen
                        )
                    ),
                    "C_OG1_A": carbon_og1,
                    "C_Nleaving_A": carbon_nleave,
                    "Nalpha_NH_bond_lengths_A": json.dumps(nalpha_nh, sort_keys=True),
                    "Nleaving_NH_bond_lengths_A": json.dumps(nleave_nh, sort_keys=True),
                }
                for name, (ha, da, dha, acceptor_id) in hbond.items():
                    row[f"waterO_{name}_D_A"] = da
                    row[f"waterH_{name}_A_A"] = ha
                    row[f"waterO_H_{name}_angle_deg"] = dha
                    row[f"waterH_{name}_closest_acceptor_atom_index_1based"] = (
                        acceptor_id
                    )
                for residue_index, atom_ids in ASP_RESIDUES.items():
                    row[f"O_Asp{residue_index}_nearestO_A"] = float(
                        np.min(
                            np.linalg.norm(
                                image_near(
                                    xyz[[atom - 1 for atom in atom_ids]], oxygen, box
                                )
                                - oxygen,
                                axis=1,
                            )
                        )
                    )
                writer.writerow(row)
    return nearest_rows, minima
