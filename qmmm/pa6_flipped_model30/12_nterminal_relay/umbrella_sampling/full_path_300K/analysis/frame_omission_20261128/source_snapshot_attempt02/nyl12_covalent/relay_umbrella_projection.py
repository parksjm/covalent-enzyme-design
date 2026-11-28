"""Equal-Angstrom path projection with explicit end and ambiguity diagnostics."""

from __future__ import annotations

from typing import Any

import numpy as np

KB = 0.00198720425864083


def project_path(points: np.ndarray, path: np.ndarray) -> dict[str, np.ndarray]:
    """Project three-CV samples onto a polyline with extended terminal tangents.

    Args:
        points: Samples in q_acyl, q_act, q_del order, in Angstrom.
        path: Ordered accepted reference images in the same units and order.

    Returns:
        Raw normalized arclength, transverse distance, nearest segment, and
        distance margin to the next nonadjacent segment. Negative margins are
        impossible; small margins indicate projection ambiguity. End projections
        can lie outside [0, 1] and are retained.
    """
    points, path = np.asarray(points, float), np.asarray(path, float)
    if (
        points.ndim != 2
        or path.ndim != 2
        or points.shape[1:] != (3,)
        or path.shape[1:] != (3,)
    ):
        raise ValueError("Samples and path must have shape (N, 3)")
    if len(path) < 2 or not np.isfinite(points).all() or not np.isfinite(path).all():
        raise ValueError("Finite samples and at least two finite path points required")
    vectors = np.diff(path, axis=0)
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths <= 0):
        raise ValueError("Consecutive reference images must differ")
    along = np.einsum("nsi,si->ns", points[:, None] - path[:-1], vectors) / lengths**2
    for segment in range(len(vectors)):
        along[:, segment] = np.clip(
            along[:, segment],
            -np.inf if segment == 0 else 0,
            np.inf if segment == len(vectors) - 1 else 1,
        )
    projected = path[:-1] + along[..., None] * vectors
    distances = np.linalg.norm(points[:, None] - projected, axis=2)
    nearest = distances.argmin(axis=1)
    rows = np.arange(len(points))
    arclength = np.r_[0, np.cumsum(lengths)]
    s = (arclength[nearest] + along[rows, nearest] * lengths[nearest]) / arclength[-1]
    nonadjacent = np.abs(np.arange(len(vectors))[None] - nearest[:, None]) > 1
    alternative = np.min(np.where(nonadjacent, distances, np.inf), axis=1)
    return dict(
        s=s,
        transverse_A=distances[rows, nearest],
        segment=nearest,
        ambiguity_margin_A=alternative - distances[rows, nearest],
    )


def path_histogram(
    s: np.ndarray,
    weights: np.ndarray,
    edges: np.ndarray,
    temperature: float,
    es_bounds: tuple[float, float],
    acyl_bounds: tuple[float, float],
) -> dict[str, Any]:
    """Estimate an unbiased histogram PMF, preserving empty bins and ES reference.

    Args:
        s: Raw path projections, including terminal fluctuations.
        weights: Common-normalization unbiased MBAR probabilities.
        edges: Explicit histogram edges covering all projections.
        temperature: Physical temperature in kelvin.
        es_bounds: Declared interval for the E-S minimum search.
        acyl_bounds: Declared interval for the acylenzyme minimum search.

    Returns:
        Density, PMF, raw per-bin sample counts, support and the contiguous
        E-S-to-acylenzyme maximum when resolved.
    """
    s = np.asarray(s, float)
    weights = np.asarray(weights, float)
    edges = np.asarray(edges, float)
    if (
        s.ndim != 1
        or weights.ndim != 1
        or len(s) == 0
        or len(weights) != len(s)
        or not np.isfinite(s).all()
        or not np.isfinite(weights).all()
        or np.any(weights < 0)
        or not np.isfinite(weights.sum())
        or weights.sum() <= 0
    ):
        raise ValueError(
            "Finite one-dimensional samples and nonnegative weights required"
        )
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be finite and positive")
    if (
        edges.ndim != 1
        or len(edges) < 3
        or not np.isfinite(edges).all()
        or np.any(np.diff(edges) <= 0)
        or s.min() < edges[0]
        or s.max() > edges[-1]
    ):
        raise ValueError("Increasing finite edges must cover all raw s support")
    if not es_bounds[0] < es_bounds[1] or not acyl_bounds[0] < acyl_bounds[1]:
        raise ValueError("Declare a nonempty E-S search interval")
    centers = (edges[:-1] + edges[1:]) / 2
    bin_indices = np.searchsorted(edges, s, side="right") - 1
    bin_indices[s == edges[-1]] = len(edges) - 2
    raw_sample_count = np.bincount(bin_indices, minlength=len(edges) - 1)
    mass = np.bincount(bin_indices, weights=weights, minlength=len(edges) - 1)
    density = mass / np.diff(edges)
    occupied = density > 0
    pmf = np.full(len(mass), np.nan)
    pmf[occupied] = -KB * temperature * np.log(density[occupied])
    reference = occupied & (centers >= es_bounds[0]) & (centers <= es_bounds[1])
    if not reference.any():
        return dict(
            status="no_ES_reference_support",
            edges=edges.tolist(),
            probability_mass=mass.tolist(),
            raw_sample_count=raw_sample_count.tolist(),
            pmf_kcal_mol=None,
        )
    zero_index = np.flatnonzero(reference)[np.argmin(pmf[reference])]
    pmf[occupied] -= pmf[zero_index]
    acyl = occupied & (centers >= acyl_bounds[0]) & (centers <= acyl_bounds[1])
    acyl_index = np.flatnonzero(acyl)[np.argmin(pmf[acyl])] if acyl.any() else None
    es_candidates = np.flatnonzero(
        (centers >= es_bounds[0]) & (centers <= es_bounds[1])
    )
    acyl_candidates = np.flatnonzero(
        (centers >= acyl_bounds[0]) & (centers <= acyl_bounds[1])
    )
    reference_resolved = (
        len(es_candidates) >= 3
        and len(acyl_candidates) >= 3
        and zero_index not in (es_candidates[0], es_candidates[-1])
        and acyl_index not in (acyl_candidates[0], acyl_candidates[-1])
        and occupied[max(0, zero_index - 1) : zero_index + 2].all()
        and acyl_index is not None
        and occupied[max(0, acyl_index - 1) : acyl_index + 2].all()
    )
    connected = (
        reference_resolved
        and acyl_index is not None
        and acyl_index > zero_index
        and occupied[zero_index : acyl_index + 1].all()
    )
    barrier_index = (
        zero_index + int(np.argmax(pmf[zero_index : acyl_index + 1]))
        if connected
        else None
    )
    return dict(
        status=(
            "supported"
            if connected
            else "reference_unresolved_no_global_barrier"
            if not reference_resolved
            else "gapped_no_global_barrier"
        ),
        edges=edges.tolist(),
        centers=centers.tolist(),
        probability_mass=mass.tolist(),
        raw_sample_count=raw_sample_count.tolist(),
        density=density.tolist(),
        occupied=occupied.tolist(),
        pmf_kcal_mol=[float(v) if np.isfinite(v) else None for v in pmf],
        reference_s=float(centers[zero_index]),
        es_search_bounds=list(es_bounds),
        acyl_search_bounds=list(acyl_bounds),
        reference_resolved=bool(reference_resolved),
        acyl_reference_s=float(centers[acyl_index]) if acyl_index is not None else None,
        basin_minimum_difference_kcal_mol=float(pmf[acyl_index])
        if acyl_index is not None
        else None,
        barrier_s=float(centers[barrier_index]) if connected else None,
        minimum_to_barrier_kcal_mol=float(pmf[zero_index : acyl_index + 1].max())
        if connected
        else None,
    )
