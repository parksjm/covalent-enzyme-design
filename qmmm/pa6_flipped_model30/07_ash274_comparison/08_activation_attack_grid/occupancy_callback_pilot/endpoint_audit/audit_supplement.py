"""Add amide-N pyramidalization and explicit donor-H geometry from restarts."""

from __future__ import annotations

import json
import logging
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

from audit_endpoints import (
    BASE,
    HASHES,
    OUT,
    SOURCE,
    angle,
    coordinates,
    distance,
    load_model,
    read_json,
    track,
)


def planar_center(x: np.ndarray, center: int, neighbors: list[int]) -> dict[str, Any]:
    """Return explicit three-neighbor center pyramidalization measures."""
    vectors = x[np.array(neighbors) - 1] - x[center - 1]
    normal = np.cross(vectors[1] - vectors[0], vectors[2] - vectors[0])
    return {
        "center_1based": center,
        "neighbors_1based": neighbors,
        "signed_volume_A3": float(np.linalg.det(vectors) / 6),
        "signed_height_A": float(-vectors[0] @ normal / np.linalg.norm(normal)),
        "neighbor_angles_deg": {
            f"{a}-{b}": angle(x, a, center, b) for a, b in combinations(neighbors, 2)
        },
        "angle_sum_deg": sum(
            angle(x, a, center, b) for a, b in combinations(neighbors, 2)
        ),
    }


def main() -> None:
    """Write supplementary geometry and matched energy reference comparisons."""
    model = load_model()
    audit = read_json(OUT / "audit.json")
    track(Path(__file__))
    neighbors = sorted(a.idx + 1 for a in model.topology.atoms[19932].bond_partners)
    assert neighbors == [19920, 19934, 19941]
    first_a = [
        json.loads(line)
        for line in track(BASE / "pilot/point_A/accepted_iterations.jsonl")
        .read_text()
        .splitlines()
    ][0]
    final_a = audit["points"]["A"]["energy"]
    result: dict[str, Any] = {
        "scope": "Retained-data supplementary audit; no engine calls",
        "points": {},
    }
    for name, record in audit["points"].items():
        path = (
            SOURCE / "optimized.ncrst"
            if name == "original_ES"
            else Path(record["endpoint_path"])
        )
        x = coordinates(path)
        hbonds = {}
        for donor, hydrogens, acceptor in [
            (1638, [1639], 19921),
            (2732, [2733, 2734], 19921),
            (13905, [13906, 13907, 13917], 20317),
            (20317, [20318, 20319], 13916),
            (14579, [14580], 14560),
        ]:
            for hydrogen in hydrogens:
                key = f"{donor}-{hydrogen}-{acceptor}"
                hbonds[key] = {
                    "H_acceptor_A": distance(x, hydrogen, acceptor),
                    "donor_acceptor_A": distance(x, donor, acceptor),
                    "donor_H_acceptor_angle_deg": angle(x, donor, hydrogen, acceptor),
                }
        point = {
            "amide_N": planar_center(x, 19933, neighbors),
            "donor_H_acceptor_geometry": hbonds,
        }
        if name != "original_ES":
            point["energy_references_kcal_mol"] = {}
            for label, ref in [("first_A", first_a), ("converged_A", final_a)]:
                for physical, key in [
                    ("physical", "physical_energy_kcal_mol"),
                    ("biased", "total_energy_kcal_mol"),
                ]:
                    point["energy_references_kcal_mol"][
                        f"{physical}_relative_{label}"
                    ] = record["energy"][key] - ref[key]
        result["points"][name] = point
    result["source_sha256"] = HASHES
    with (OUT / "supplement.json").open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    logging.info(
        "Supplemented amide-N pyramidalization, donor-H angles, biased A references"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
