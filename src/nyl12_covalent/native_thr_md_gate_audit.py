"""Diagnostic analysis of the original four native-Thr geometry gates."""

from __future__ import annotations

from collections import Counter
from typing import Any
import numpy as np

KEYS = ("og1_shift", "nterm_shift", "psi_delta", "chi1_delta")
LIMITS = np.array([0.75, 1.25, 35.0, 20.0])
GRID = (
    (0.75, 1.0, 1.25, 1.5),
    (1.25, 1.5, 2.0, 2.5),
    (35.0, 45.0, 60.0, 90.0),
    (20.0, 30.0, 40.0, 60.0),
)


def values(rows: list[dict[str, Any]]) -> np.ndarray:
    """Extract finite absolute gate metrics without changing original limits.

    Args:
        rows: Saved geometry records.

    Returns:
        Array of four absolute gate metrics per row.
    """
    matrix = np.array([[abs(float(r[k])) for k in KEYS] for r in rows], dtype=float)
    if matrix.shape != (len(rows), 4) or not np.isfinite(matrix).all():
        raise ValueError("Invalid gate metrics")
    return matrix


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize gate failures, combinations and diagnostic sensitivity.

    Args:
        rows: Nonempty saved frame/site records.

    Returns:
        Counts and distributions with an explicit denominator.
    """
    if not rows:
        raise ValueError("No rows")
    x = values(rows)
    fail = x > LIMITS
    counts = Counter(
        "+".join(KEYS[i] for i in range(4) if flags[i]) or "pass" for flags in fail
    )
    stats = {}
    for key in (*KEYS, "scaffold_rmsd"):
        a = np.array(
            [
                abs(float(r[key]))
                if key in ("psi_delta", "chi1_delta")
                else float(r[key])
                for r in rows
            ]
        )
        if not np.isfinite(a).all():
            raise ValueError("Invalid scaffold metric")
        stats[key] = {
            name: float(val)
            for name, val in zip(
                ("q05", "median", "q95"), np.quantile(a, [0.05, 0.5, 0.95]), strict=True
            )
        }
    single = {}
    for i, key in enumerate(KEYS):
        single[key] = [
            {
                "limit": float(limit),
                "passes": int(
                    np.all(
                        x <= np.where(np.arange(4) == i, limit, LIMITS), axis=1
                    ).sum()
                ),
            }
            for limit in GRID[i]
        ]
    joint = []
    for bounds in (
        (1.0, 1.5, 45.0, 30.0),
        (1.25, 2.0, 60.0, 40.0),
        (1.5, 2.5, 90.0, 60.0),
    ):
        joint.append(
            {
                "limits": dict(zip(KEYS, bounds, strict=True)),
                "passes": int(np.all(x <= bounds, axis=1).sum()),
            }
        )
    return {
        "n": len(rows),
        "fails": dict(zip(KEYS, fail.sum(axis=0).tolist(), strict=True)),
        "fail_percent": dict(
            zip(KEYS, (100 * fail.mean(axis=0)).tolist(), strict=True)
        ),
        "failure_combinations": dict(counts),
        "exclusive_fails": {
            key: int(np.sum(fail[:, i] & (fail.sum(axis=1) == 1)))
            for i, key in enumerate(KEYS)
        },
        "leave_one_out_passes": {
            key: int(np.all(~np.delete(fail, i, axis=1), axis=1).sum())
            for i, key in enumerate(KEYS)
        },
        "distributions": stats,
        "one_at_a_time_diagnostic": single,
        "joint_diagnostic": joint,
    }


def representatives(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Choose deterministic Tier-A/B AE median medoid and closest-native frame.

    Args:
        rows: Saved records from a single replica.

    Returns:
        Named source rows selected by predefined normalized-distance rules.
    """
    eligible = [r for r in rows if r["tier"] in ("A", "B")]
    ae = [r for r in eligible if r["site"] == "AE"]
    if not ae or not eligible:
        raise ValueError("No Tier A/B frames for selection")
    x = values(ae) / LIMITS
    median = np.median(x, axis=0)
    medoid = min(
        range(len(ae)),
        key=lambda i: (float(np.linalg.norm(x[i] - median)), float(ae[i]["time_ps"])),
    )
    scores = values(eligible) / LIMITS
    closest = min(
        range(len(eligible)),
        key=lambda i: (
            float(np.max(scores[i])),
            float(eligible[i]["time_ps"]),
            eligible[i]["site"],
        ),
    )
    return {"AE_median_medoid": ae[medoid], "closest_native": eligible[closest]}
