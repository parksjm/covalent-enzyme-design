"""Smooth algebraic internal-coordinate path and stable local histograms."""

from __future__ import annotations

import numpy as np
from scipy.special import logsumexp


def uniform_nodes(path: np.ndarray, count: int = 153) -> tuple[np.ndarray, float]:
    """Resample a nondegenerate polyline at fixed uniform arclength positions."""
    lengths = np.linalg.norm(np.diff(path, axis=0), axis=1)
    if np.any(lengths <= 0) or count < 2:
        raise ValueError("Degenerate path or node count")
    arc = np.r_[0.0, np.cumsum(lengths)]
    target = np.linspace(0, arc[-1], count)
    return np.column_stack([np.interp(target, arc, x) for x in path.T]), arc[-1]


def smooth_progress(
    q: np.ndarray, nodes: np.ndarray, lam: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return BGP normalized progress and its analytic gradient in q space."""
    if lam <= 0 or not np.isfinite(q).all():
        raise ValueError("Invalid lambda or coordinates")
    u = np.linspace(0, 1, len(nodes))
    delta = q[:, None, :] - nodes[None, :, :]
    logits = -lam * np.sum(delta * delta, axis=2)
    probability = np.exp(logits - logsumexp(logits, axis=1)[:, None])
    s = probability @ u
    gradient = (
        -2
        * lam
        * np.einsum("ij,ij,ijk->ik", probability, u[None, :] - s[:, None], delta)
    )
    return s, gradient


def local_histogram(
    s: np.ndarray, logw: np.ndarray, edges: np.ndarray
) -> dict[str, list]:
    """Accumulate each bin in log space and preserve unsupported bins as null."""
    index = np.searchsorted(edges, s, side="right") - 1
    index[s == edges[-1]] = len(edges) - 2
    normalized = logw - logsumexp(logw)
    density, mass, ess, counts = [], [], [], []
    for i, width in enumerate(np.diff(edges)):
        local = normalized[index == i]
        counts.append(int(len(local)))
        if len(local):
            lm = logsumexp(local)
            density.append(float(lm - np.log(width)))
            mass.append(float(np.exp(lm)))
            ess.append(float(np.exp(2 * lm - logsumexp(2 * local))))
        else:
            density.append(None)
            mass.append(0.0)
            ess.append(None)
    return dict(logdensity=density, mass=mass, ess=ess, counts=counts)
