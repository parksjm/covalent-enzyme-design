"""Separate heavy-atom barriers from hydrogen barriers; no acceptance widening."""

from __future__ import annotations

import json
import logging

import numpy as np
from search_oxygen import OUT, front_frame, load_pocket, screen_oxygen


def contact_table(point: np.ndarray, environment: tuple) -> list[dict]:
    """Return the closest environment contacts with explicit vdW margins.

    Args:
        point: Diagnostic oxygen coordinate.
        environment: Complete fixed pocket environment.

    Returns:
        Twenty contacts ranked by unexceptioned vdW overlap.
    """
    xyz, radius, names, _ = environment
    distance = np.linalg.norm(xyz - point, axis=1)
    margin = distance - (1.52 + radius - 0.4)
    rank = np.argsort(margin)[:20]
    return [
        {
            "atom": names[i],
            "distance_A": float(distance[i]),
            "unexceptioned_margin_A": float(margin[i]),
        }
        for i in rank
    ]


def main() -> None:
    """Write diagnostic heavy-only and widened-distance oxygen results."""
    _, seed, _, env = load_pocket()
    grid = np.load(OUT / "oxygen_grid.npz")
    points, front = grid["coordinates"], grid["front"]
    heavy = np.flatnonzero(env[1] != 1.2)
    heavy_env = (
        env[0][heavy],
        env[1][heavy],
        [env[2][i] for i in heavy],
        env[3][heavy],
    )
    margin, collider, _ = screen_oxygen(points, heavy_env, seed[0])
    preferred = np.array([-3.833305, 8.999226, -1.874442])
    rank = np.lexsort((np.linalg.norm(points - preferred, axis=1), -margin))[:10]
    report = {
        "heavy_only_acceptance": False,
        "reason": "Diagnostic only; every H contact remains required for acceptance",
        "heavy_only_pass_count": int(np.sum(margin >= -1e-8)),
        "heavy_only_front_sensitivity": {
            str(c): int(np.sum((margin >= -1e-8) & (front >= c)))
            for c in [0, 0.25, 0.5, 0.75]
        },
        "heavy_best": [
            {
                "oxygen_A": points[i].tolist(),
                "heavy_margin_A": float(margin[i]),
                "heavy_collider": heavy_env[2][collider[i]],
                "front_projection": float(front[i]),
                "full_contacts": contact_table(points[i], env),
            }
            for i in rank
        ],
        "accepted_range_best_allatom_contacts": contact_table(
            points[np.argmax(grid["margin"])], env
        ),
        "user_oxygen_contacts": contact_table(preferred, env),
    }
    # Independently declared wider grid: not an accepted replacement range.
    axis, forward, right = front_frame(seed[0], seed[3], preferred)
    distance = np.linalg.norm(seed[3] - seed[0])
    rn, ro, phi = np.meshgrid(
        np.linspace(2.6, 3.3, 29),
        np.linspace(2.6, 3.3, 29),
        np.deg2rad(np.arange(-89.75, 90, 0.5)),
        indexing="ij",
    )
    x = (rn**2 - ro**2 + distance**2) / (2 * distance)
    rho = np.sqrt(rn**2 - x**2)
    wide = (
        seed[0]
        + x[..., None] * axis
        + rho[..., None]
        * (np.cos(phi)[..., None] * forward + np.sin(phi)[..., None] * right)
    ).reshape(-1, 3)
    wide_margin, wide_collider, _ = screen_oxygen(
        wide,
        env,
        seed[0],
        max_n_distance=3.3,
    )
    best = np.argmax(wide_margin)
    report["widened_3p3_diagnostic"] = {
        "accepted": False,
        "point_count": len(wide),
        "pass_count": int(np.sum(wide_margin >= -1e-8)),
        "best_oxygen_A": wide[best].tolist(),
        "best_margin_A": float(wide_margin[best]),
        "best_collider": env[2][wide_collider[best]],
        "contacts": contact_table(wide[best], env),
    }
    (OUT / "oxygen_diagnostics.json").write_text(json.dumps(report, indent=2) + "\n")
    logging.warning(
        "Heavy-only accepted-grid passes %d; widened allatom passes %d",
        report["heavy_only_pass_count"],
        np.sum(wide_margin >= -1e-8),
    )


if __name__ == "__main__":
    main()
