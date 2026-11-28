"""Replica-level occupancy statistics for the matched twelve-run WT cohort."""

from __future__ import annotations

from collections import Counter
import math
from typing import Any

import numpy as np

from nyl12_covalent.modular_es_snapshot_selection import classify_geometry

STATES = {"nh2_asp272_asp274", "nh2_asp272_ash274", "nh3_asp272_asp274"}
OBSERVABLES = (
    "tier_a",
    "tier_b_exclusive",
    "tier_a_or_b",
    "native",
    "tier_a_and_native",
    "tier_b_exclusive_and_native",
    "tier_a_or_b_and_native",
    "phi_positive",
    "phi_negative",
    "phi_zero",
    "tier_a_phi_positive",
    "tier_a_phi_negative",
    "tier_a_or_b_phi_positive",
    "tier_a_or_b_phi_negative",
)
NATIVE_GATES = {
    "og1_shift_A": 0.75,
    "nalpha_shift_A": 1.25,
    "psi_delta_deg": 35.0,
    "chi1_delta_deg": 20.0,
}


def validate_config(config: dict[str, Any]) -> None:
    """Require the matched three-state, four-independent-replica manifest.

    Args:
        config: JSON configuration with explicit cases, cadence and native gates.
    """
    cases = config["cases"]
    if len(cases) != 12 or Counter(c["state"] for c in cases) != Counter(
        {state: 4 for state in STATES}
    ):
        raise ValueError("Require exactly the three matched WT states, four runs each")
    names = [c["name"] for c in cases]
    if len(set(names)) != 12 or any(
        not name
        or any(
            c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
            for c in name
        )
        for name in names
    ):
        raise ValueError("Require twelve unique safe case names")
    if float(config["dt_ps"]) != 10 or config["native_gates"] != NATIVE_GATES:
        raise ValueError("Require matched 10 ps cadence and established native gates")
    if float(config["time_offset_ps"]) != 500:
        raise ValueError("Require the declared 500 ps preparation offset")
    paths = [path for case in cases for path in case["trajectories"]]
    if len(set(paths)) != len(paths):
        raise ValueError("Trajectory paths must be independent across cases")
    for case in cases:
        if case["expected_frames"] != 10000 or not case["trajectories"]:
            raise ValueError("Require all 10000 production frames per case")
        if not {
            "og1",
            "nalpha",
            "carbon",
            "oxygen",
            "nleave",
            "adjacent_c",
            "oxyanion1",
            "oxyanion2",
        }.issubset(case["geometry_atoms"]):
            raise ValueError("Explicit complete case-owned geometry mapping required")


def frame_flags(metrics: dict[str, float], native: bool) -> dict[str, bool]:
    """Classify exclusive tiers, native intersections and independent phi signs.

    Args:
        metrics: Fresh finite catalytic geometry.
        native: Whether all held-out native Thr gates pass.

    Returns:
        Boolean observables counted using all production frames as denominator.
    """
    if not all(math.isfinite(value) for value in metrics.values()):
        raise ValueError("Nonfinite catalytic geometry")
    tier = classify_geometry(metrics)
    a, b = tier == "A", tier == "B"
    positive, negative = metrics["projected_phi"] > 0, metrics["projected_phi"] < 0
    return dict(
        zip(
            OBSERVABLES,
            (
                a,
                b,
                a or b,
                native,
                a and native,
                b and native,
                (a or b) and native,
                positive,
                negative,
                not (positive or negative),
                a and positive,
                a and negative,
                (a or b) and positive,
                (a or b) and negative,
            ),
            strict=True,
        )
    )


def block_index(frame: int, dt_ps: float) -> int:
    """Return one-based 10 ns block for saved production endpoints.

    Args:
        frame: One-based ordinal production frame (preparation offset excluded).
        dt_ps: Production sampling interval in ps.

    Returns:
        Block with right-closed production endpoints, e.g. 10 ns in block 1.
    """
    if frame < 1 or dt_ps <= 0:
        raise ValueError("Require positive frame and cadence")
    return int(math.ceil(frame * dt_ps / 10000))


def fractions(counts: Counter[str], frames: int) -> dict[str, float]:
    """Normalize observable counts by all saved frames.

    Args:
        counts: Integer observable counts.
        frames: Number of frames in this run or block.

    Returns:
        Explicit zero-inclusive occupancy fractions.
    """
    if frames <= 0:
        raise ValueError("Cannot summarize empty frames")
    return {key: counts[key] / frames for key in OBSERVABLES}


def summarize_states(replicas: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize equally weighted independent replicas, without frame-level CIs.

    Args:
        replicas: Four per-state summaries with complete-frame occupancies.

    Returns:
        Mean, sample SD and range for each observable, with n=4 explicitly.
    """
    if Counter(r["state"] for r in replicas) != Counter({s: 4 for s in STATES}):
        raise ValueError("State summaries require four independent replicas each")
    return {
        state: {
            "replicas": 4,
            "observables": {
                key: {
                    "mean": float(np.mean(values)),
                    "sample_sd": float(np.std(values, ddof=1)),
                    "min": float(min(values)),
                    "max": float(max(values)),
                }
                for key in OBSERVABLES
                for values in [
                    [r["fractions"][key] for r in replicas if r["state"] == state]
                ]
            },
        }
        for state in sorted(STATES)
    }
