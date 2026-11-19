"""Reproduce recent-band peak/neighborhood extraction from verified retained inputs."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


def extract(manifest: Path) -> dict[str, Any]:
    """Verify frozen source hashes and extract peak neighbors and H-owner changes.

    Args:
        manifest: Retained diagnosis JSON declaring exact source file SHA-256s.

    Returns:
        Recomputed scalar geometry, proton inventory and owner-change evidence.
    """
    source = json.loads(manifest.read_text())
    mapping = source["atom_map_1based"]
    peaks, neighbors, changes, previous = [], [], [], {}
    for name, expected in source["input_sha256"].items():
        path = Path(name)
        payload = path.read_bytes()
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError(f"Retained diagnosis source changed: {path}")
        if path.name != "images.json":
            continue
        images = json.loads(payload)
        peak = max(images, key=lambda row: row["relative_energy_kcal_mol"])
        cycle = peak["cycle"]
        peaks.append(
            {
                "cycle": cycle,
                "image": peak["image"],
                "maximum_kcal_mol": peak["relative_energy_kcal_mol"],
            }
        )
        for row in images:
            owners = {h["hydrogen_1based"]: h for h in row["state"]["proton_owners"]}
            image = row["image"]
            if image in previous:
                old_cycle, old = previous[image]
                moved = [
                    {
                        "H": h,
                        "from": old[h]["nearest_qm_heavy_1based"],
                        "to": value["nearest_qm_heavy_1based"],
                    }
                    for h, value in owners.items()
                    if old[h]["nearest_qm_heavy_1based"]
                    != value["nearest_qm_heavy_1based"]
                ]
                if moved:
                    changes.append(
                        {
                            "from_cycle": old_cycle,
                            "cycle": cycle,
                            "image": image,
                            "owner_changes": moved,
                        }
                    )
            previous[image] = cycle, owners
            if abs(image - peak["image"]) <= 1:
                counts = Counter(
                    h["nearest_qm_heavy_1based"]
                    for h in owners.values()
                    if h["nearest_distance_A"] <= 1.25
                )
                neighbors.append(
                    {
                        "cycle": cycle,
                        "image": image,
                        "relative_energy_kcal_mol": row["relative_energy_kcal_mol"],
                        "geometry": {
                            k: v
                            for k, v in row["state"].items()
                            if isinstance(v, (int, float))
                        },
                        "counts": {
                            k: counts[mapping[k]]
                            for k in (
                                "thr_n",
                                "thr_og",
                                "water_o",
                                "leaving_n",
                                "asp_o1",
                                "asp_o2",
                                "ash_oxygen",
                            )
                        },
                        "special_protons": {
                            k: owners[mapping[k]]
                            for k in ("proton", "water_h1", "water_h2", "ash_proton")
                        },
                    }
                )
    return {
        "maxima": peaks,
        "maxima_and_neighbors": neighbors,
        "all_image_nearest_H_owner_changes": changes,
    }


if __name__ == "__main__":
    path = Path(__file__).with_name("recent_band_geometry_diagnosis_20261119.json")
    json.dump(extract(path), sys.stdout, indent=2)
    sys.stdout.write("\n")
