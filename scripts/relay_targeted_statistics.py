"""Small, explicit statistical screens for targeted relay sampling."""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from pymbar import timeseries

FEATURES = ("q_acyl_A", "q_act_A", "q_del_A", "r_c_og1", "r_c_n")
OWNERS = ("hg1_owner", "nalpha_h1_owner", "nalpha_h2_owner")


def production_rows(
    paths: list[Path], discard_ps: float = 5.0
) -> tuple[list[dict], list[dict]]:
    """Concatenate ordered production blocks and discard once per trajectory.

    Args:
        paths: Chronologically ordered production CSV paths (equilibration excluded).
        discard_ps: Fixed initial production interval to remove, exactly once.

    Returns:
        All frames and retained frames with elapsed production times.
    """
    rows = []
    offset_steps = 0
    origin_steps = None
    for path in paths:
        mdin = path.with_name("run.mdin").read_text()

        def parameter(name: str) -> float:
            match = re.search(r"\b" + name + r"\s*=\s*([-+\d.eE]+)", mdin, re.I)
            if match is None:
                raise ValueError(f"Missing declared {name}: {path}")
            return float(match[1])

        dt = parameter("dt")
        stride, total = int(parameter("ntwx")), int(parameter("nstlim"))
        if dt != 0.0005 or stride <= 0 or total % stride:
            raise ValueError("Unexpected declared integration clock or cadence")
        with path.open() as handle:
            block = list(csv.DictReader(handle))
        if len(block) < 2 or len(block) != total // stride:
            raise ValueError(
                "Require at least two frames and complete declared cadence"
            )
        for index, row in enumerate(block, 1):
            if int(row["frame"]) != index:
                raise ValueError("Noncontiguous local frame indices")
            local_step = index * stride
            native = float(row["time_ps"])
            if origin_steps is None:
                origin_steps = round(native / dt) - local_step
            expected = (origin_steps + offset_steps + local_step) * dt
            tolerance = float(abs(np.spacing(np.float32(expected))))
            if abs(native - expected) > tolerance + 1e-10:
                raise ValueError("Native clock discontinuity beyond float32 ULP")
            rows.append(dict(row, production_time_ps=(offset_steps + local_step) * dt))
        offset_steps += total
    if not rows:
        raise ValueError("No production frames")
    return rows, [r for r in rows if r["production_time_ps"] > discard_ps]


def summarize(rows: list[dict]) -> dict[str, Any]:
    """Return uncensored feature, correlation, drift and proton-owner screens.

    Args:
        rows: Retained contiguous trajectory frames.

    Returns:
        Counts, feature statistics and complete tracked-proton owner frequencies.
    """
    if len(rows) < 4:
        raise ValueError("At least four retained frames required")
    result: dict[str, Any] = {"frames": len(rows), "features": {}, "owners": {}}
    for name in FEATURES:
        x = np.array([float(r[name]) for r in rows])
        if not np.all(np.isfinite(x)):
            raise ValueError("Nonfinite scientific feature")
        constant = bool(np.ptp(x) < 1e-12)
        g = float(timeseries.statistical_inefficiency(x)) if not constant else None
        half = len(x) // 2
        result["features"][name] = dict(
            mean=float(x.mean()),
            sd=float(x.std(ddof=1)),
            minimum=float(x.min()),
            maximum=float(x.max()),
            half_drift_A=float(x[half:].mean() - x[:half].mean()),
            statistical_inefficiency=g,
            ess=None if g is None else len(x) / g,
            constant=constant,
            constant_note="ESS unresolved for constant observable"
            if constant
            else None,
        )
    for name in OWNERS:
        result["owners"][name] = dict(Counter(r[name] for r in rows))
    result["hydrogen_inventory"] = dict(Counter(r["hydrogen_owners"] for r in rows))
    return result


def histogram_overlap(a: np.ndarray, b: np.ndarray, width: float = 0.1) -> float:
    """Calculate histogram intersection on pooled, zero-anchored 1D/2D bins.

    Args:
        a: First samples, optionally two dimensional.
        b: Second samples of matching dimension.
        width: Bin width in Angstrom.

    Returns:
        Probability histogram intersection between zero and one.
    """
    a, b = np.asarray(a), np.asarray(b)
    if a.ndim == 1:
        a, b = a[:, None], b[:, None]
    pooled = np.concatenate((a, b))
    bins = [
        np.arange(
            np.floor(x.min() / width) * width,
            (np.ceil(x.max() / width) + 2) * width,
            width,
        )
        for x in pooled.T
    ]
    ha = np.histogramdd(a, bins=bins)[0] / len(a)
    hb = np.histogramdd(b, bins=bins)[0] / len(b)
    return float(np.minimum(ha, hb).sum())


def require_same_bias(a: dict, b: dict) -> None:
    """Reject same-bias comparisons unless full native bias definitions match.

    Args:
        a: First window bias definition.
        b: Second window bias definition.
    """
    if a != b:
        raise ValueError("Reciprocal comparison requires identical native bias")


def recommendation(
    ess: float | None,
    drift: float,
    overlap: float,
    same_bias: bool,
    biased_overlap: float,
    gap: float,
    chemistry_failure: bool = False,
) -> str:
    """Apply prospective diagnostic heuristics without claiming convergence.

    Args:
        ess: Minimum resolved feature effective sample count.
        drift: Largest absolute half-trajectory drift in Angstrom.
        overlap: Physical q_del histogram intersection.
        same_bias: Whether native biases are identical.
        biased_overlap: Joint biased-CV histogram intersection.
        gap: Absolute q_del mean gap in Angstrom.
        chemistry_failure: Explicit reviewed chemistry failure.

    Returns:
        Conditional next-step recommendation; never a launch instruction.
    """
    if chemistry_failure:
        return "inspect_chemistry_before_more_sampling"
    if ess is None or ess < 20 or drift > 0.2:
        return "longer_sampling_required_before_center_decision"
    if same_bias and biased_overlap >= 0.5 and (overlap < 0.5 or gap > 0.3):
        return "orthogonal_qdel_bridge_or_multidimensional_sampling_candidate"
    if overlap < 0.2:
        return "additional_centers_candidate_if_poor_overlap_persists"
    return "conditional_no_new_center_evidence_extend_product_toward_30ps"
