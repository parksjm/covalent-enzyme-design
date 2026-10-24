"""Opposite-hemisphere oxygen diagnostic using the unchanged pocket screens."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
from diagnose_oxygen import contact_table
from search_oxygen import OUT, front_frame, load_pocket, screen_oxygen


def reflect_radial(points: np.ndarray, n: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Rotate radial projections by 180 degrees about the N-OG axis.

    Args:
        points: Front-side oxygen coordinates.
        n: Nitrogen coordinate.
        axis: Unit N-to-OG vector.

    Returns:
        Opposite-side points with identical distances to N and OG.
    """
    displacement = points - n
    axial = np.outer(displacement @ axis, axis)
    return n + 2 * axial - displacement


def main() -> None:
    """Preserve opposite oxygen grid and steric necessary-condition results."""
    _, seed, _, environment = load_pocket()
    grid_path = OUT / "oxygen_grid.npz"
    original = np.load(grid_path)
    preferred = np.array([-3.833305, 8.999226, -1.874442])
    axis, _, _ = front_frame(seed[0], seed[3], preferred)
    points = reflect_radial(original["coordinates"], seed[0], axis)
    front = -original["front"]
    full_margin, full_collider, _ = screen_oxygen(points, environment, seed[0])
    heavy_ids = np.flatnonzero(environment[1] != 1.2)
    heavy_env = (
        environment[0][heavy_ids],
        environment[1][heavy_ids],
        [environment[2][i] for i in heavy_ids],
        environment[3][heavy_ids],
    )
    heavy_margin, heavy_collider, _ = screen_oxygen(points, heavy_env, seed[0])
    passed = full_margin >= -1e-8
    # Rank genuinely clear oxygen candidates first by distance from user's O.
    clear_ids = np.flatnonzero(passed)
    clear_rank = clear_ids[
        np.argsort(np.linalg.norm(points[clear_ids] - preferred, axis=1))
    ][:10]
    report = {
        "diagnostic_only": True,
        "candidate_seed_emitted": False,
        "not_mechanism_evidence": True,
        "grid_count": len(points),
        "N_OG_distance_range_A": [2.6, 3.2],
        "front_projection": "negative of original normalized front radial projection",
        "full_pass_count": int(passed.sum()),
        "heavy_pass_count": int(np.sum(heavy_margin >= -1e-8)),
        "best_full_margin_A": float(full_margin.max()),
        "best_full_collider": environment[2][full_collider[np.argmax(full_margin)]],
        "best_heavy_margin_A": float(heavy_margin.max()),
        "best_heavy_collider": heavy_env[2][heavy_collider[np.argmax(heavy_margin)]],
        "full_clear_front_sensitivity": {
            str(c): int(np.sum(passed & (front <= -c))) for c in [0, 0.25, 0.5, 0.75]
        },
        "clear_oxygen_positions": [
            {
                "oxygen_A": points[i].tolist(),
                "front_projection": float(front[i]),
                "distance_from_user_O_A": float(np.linalg.norm(points[i] - preferred)),
                "full_margin_A": float(full_margin[i]),
                "contacts": contact_table(points[i], environment),
            }
            for i in clear_rank
        ],
        "source_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                Path(__file__),
                Path(__file__).with_name("search_oxygen.py"),
                Path(__file__).with_name("diagnose_oxygen.py"),
                grid_path,
            ]
        },
    }
    np.savez(
        OUT / "opposite_oxygen_grid.npz",
        coordinates=points,
        front=front,
        full_margin=full_margin,
        heavy_margin=heavy_margin,
    )
    (OUT / "opposite_oxygen_feasibility.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    logging.warning(
        "Opposite full clear %d/%d; heavy clear %d/%d",
        passed.sum(),
        len(points),
        np.sum(heavy_margin >= -1e-8),
        len(points),
    )


if __name__ == "__main__":
    main()
