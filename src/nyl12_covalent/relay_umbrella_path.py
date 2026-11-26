"""Densify accepted relay paths in their chart-specific two-CV spaces."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow
from nyl12_covalent.relay_umbrella_prepare import (
    PILOT_FORCE_CONSTANT,
    PILOT_OUTER_EXTENT_A,
    PilotWindow,
    digest,
    geometry_features,
    seed_candidates,
    select_seed,
)

def _source_candidate(base: Path) -> dict[str, Any]:
    """Build the shared source seed candidate from retained identity evidence."""
    identity_path = base / "model/source_identity.json"
    identity = json.loads(identity_path.read_text())
    acceptance = json.loads((base / "model/release_acceptance.json").read_text())
    receipt = {
        "native_converged": acceptance.get("native_converged"),
        "electronic_history": acceptance.get("electronic_history", {}),
        "geometry": {"distances_A": identity["initial_geometry"]["distances_A"]},
    }
    values = geometry_features(receipt)
    restart = base / "model/inputs/optimized.ncrst"
    if not restart.is_file():
        raise FileNotFoundError(restart)
    return {
        "seed_id": "source_w033",
        "window": "source_w033",
        "restart": restart,
        "evidence": identity_path,
        "sha256": _sha(restart),
        "features": values,
    }


def _sha(path: Path) -> str:
    """Return SHA-256 for a seed file."""
    return digest(path)


def path_candidates(base: Path, branch: str) -> list[dict[str, Any]]:
    """Return shared source, accepted driven windows, and released endpoint."""
    if branch not in {"01_activation_reverse", "02_acylation_forward"}:
        raise ValueError(f"Unsupported relay path branch: {branch}")
    branch_dir = base / branch
    candidates = [_source_candidate(base), *seed_candidates(branch_dir)]
    endpoint_evidence = branch_dir / "02_unrestrained_endpoint/acceptance.json"
    endpoint_acceptance = json.loads(endpoint_evidence.read_text())
    values = geometry_features(endpoint_acceptance)
    endpoint_restart = branch_dir / "02_unrestrained_endpoint/optimized.ncrst"
    if Path(endpoint_acceptance["restart"]).resolve() != endpoint_restart.resolve():
        raise ValueError("Endpoint evidence points to an unexpected restart")
    candidates.append(
        {
            "seed_id": f"{branch}_endpoint",
            "window": "unrestrained_endpoint",
            "restart": endpoint_restart,
            "evidence": endpoint_evidence,
            "sha256": _sha(endpoint_restart),
            "features": values,
        }
    )
    return candidates


def _densify(
    chart: str,
    feature_names: tuple[str, str],
    route: Sequence[Mapping[str, Any]],
    candidates: Sequence[Mapping[str, Any]],
    maximum_axis_step_A: float,
) -> list[PilotWindow]:
    """Interpolate CV centers only, seeding each from a verified path restart."""
    if len(route) < 2 or maximum_axis_step_A <= 0:
        raise ValueError("Need a multi-point path and positive per-axis spacing")
    if not math.isfinite(maximum_axis_step_A):
        raise ValueError("Maximum CV spacing must be finite")
    result: list[PilotWindow] = []
    pairs: list[
        tuple[tuple[float, float], Mapping[str, Any], float, tuple[float, float]]
    ] = []
    first_values = route[0]["features"]
    first_center = tuple(float(first_values[name]) for name in feature_names)
    pairs.append((first_center, route[0], 0.0, (0.0, 0.0)))
    for left, right in zip(route, route[1:]):
        a = np.asarray([left["features"][name] for name in feature_names], dtype=float)
        b = np.asarray([right["features"][name] for name in feature_names], dtype=float)
        delta = b - a
        pieces = max(1, int(math.ceil(float(np.max(np.abs(delta))) / maximum_axis_step_A)))
        for index in range(1, pieces + 1):
            target = a + delta * (index / pieces)
            center = (float(target[0]), float(target[1]))
            if index == pieces:
                seed, distance = right, 0.0
            else:
                seed, distance = select_seed(candidates, feature_names, center)
            offset = tuple(
                float(seed["features"][name] - value)
                for name, value in zip(feature_names, center, strict=True)
            )
            bias_energy = sum(PILOT_FORCE_CONSTANT * value**2 for value in offset)
            if bias_energy > 5.0:
                raise ValueError(
                    f"Nearest accepted seed is {bias_energy:.3f} kcal/mol from {chart} center"
                )
            pairs.append((center, seed, distance, offset))
    for index, (center, seed, distance, offset) in enumerate(pairs, start=1):
        result.append(
            PilotWindow(
                label=f"{chart}_full_{index:03d}",
                chart=chart,
                bias=RelayBiasWindow(
                    feature_names=feature_names,
                    centers_A=center,
                    force_constants=(PILOT_FORCE_CONSTANT, PILOT_FORCE_CONSTANT),
                    outer_extent_A=PILOT_OUTER_EXTENT_A,
                ),
                seed_source=seed["restart"],
                seed_evidence=seed["evidence"],
                seed_sha256=seed["sha256"],
                seed_id=seed["seed_id"],
                seed_window=seed["window"],
                seed_distance_A=distance,
                seed_offset_A=offset,
            )
        )
    return result


def full_path_windows(
    base: Path, maximum_axis_step_A: float = 0.15
) -> list[PilotWindow]:
    """Create complete E–S→relay→acylenzyme CV center schedules.

    Activation follows the reversed activation path plus the shared source in
    ``(q_acyl,q_act)``. Delivery follows the shared source, a common relay
    seed, and the accepted forward path in ``(q_acyl,q_del)``. Only CV centers
    are interpolated; coordinates always come from accepted restart files.
    """
    activation = path_candidates(base, "01_activation_reverse")
    delivery = path_candidates(base, "02_acylation_forward")
    activation_route = list(reversed(activation))
    delivery_route = delivery
    specs = (
        ("activation", ("q_acyl", "q_act"), activation_route, activation),
        ("delivery", ("q_acyl", "q_del"), delivery_route, delivery),
    )
    return [
        window
        for chart, names, route, candidates in specs
        for window in _densify(
            chart, names, route, candidates, maximum_axis_step_A
        )
    ]


def equal_arclength_reference(base: Path, images: int = 39) -> list[list[float]]:
    """Return exact accepted 3-CV path nodes and cumulative arc positions.

    The analysis projection uses these nodes in the original image order and
    equal weights for all three Angstrom-valued CVs. Tangent extension belongs
    to analysis evaluation and does not add synthetic path states here.
    """
    if images != 39:
        raise ValueError("The retained accepted path contains exactly 39 images")
    reverse = path_candidates(base, "01_activation_reverse")
    forward = path_candidates(base, "02_acylation_forward")
    route = [reverse[-1], *reversed(reverse[1:-1]), reverse[0], *forward[1:]]
    points = np.asarray(
        [
            [candidate["features"][name] for name in ("q_acyl", "q_act", "q_del")]
            for candidate in route
        ],
        dtype=float,
    )
    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if not np.isfinite(segment_lengths).all() or np.any(segment_lengths <= 0):
        raise ValueError("3-CV reference path has a zero or invalid segment")
    if points.shape != (images, 3) or not np.isfinite(points).all():
        raise ValueError("Invalid equal-arclength 3-CV reference")
    return points.tolist()


def cumulative_arclength_A(base: Path) -> list[float]:
    """Return cumulative Euclidean arc length for the exact 39 path nodes."""
    points = np.asarray(equal_arclength_reference(base), dtype=float)
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))].tolist()


def maximum_center_steps(windows: Sequence[PilotWindow]) -> dict[str, float]:
    """Measure largest absolute adjacent target step on each chart axis."""
    result = {}
    for chart in ("activation", "delivery"):
        selected = [window for window in windows if window.chart == chart and "common" not in window.label]
        steps = np.diff(np.asarray([window.bias.centers_A for window in selected]), axis=0)
        result[chart] = float(np.abs(steps).max()) if len(steps) else 0.0
    return result
