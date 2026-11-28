"""Curiosity-only observation deletion and exact-bias MBAR refit helpers."""

from __future__ import annotations

import numpy as np
from pymbar import MBAR
from scipy.sparse.csgraph import connected_components
from scipy.special import logsumexp


def omit_observation(
    q: np.ndarray, counts: np.ndarray, state: int, target: np.ndarray
) -> tuple[np.ndarray, np.ndarray, int]:
    """Match one observation within its provenance state, then decrement that state.

    Args:
        q: Original state-ordered physical CV observations.
        counts: Original observations per state.
        state: Provenance-verified state index.
        target: Independently verified native-source CV coordinates.

    Returns:
        Retention mask, decremented state counts and matched original index.
    """
    counts = np.asarray(counts, int)
    if counts.ndim != 1 or counts.sum() != len(q) or np.any(counts <= 0):
        raise ValueError("Invalid original state counts")
    if not 0 <= state < len(counts) or counts[state] <= 1:
        raise ValueError("Deletion must leave a sampled source state")
    owner = np.repeat(np.arange(len(counts)), counts)
    matches = np.flatnonzero(
        (owner == state) & np.all(np.isclose(q, target, atol=1e-14, rtol=0), axis=1)
    )
    if len(matches) != 1:
        raise ValueError("Native provenance must identify exactly one observation")
    index = int(matches[0])
    keep = np.arange(len(q)) != index
    new_counts = counts.copy()
    new_counts[state] -= 1
    if new_counts.sum() != keep.sum() or keep.sum() != len(q) - 1:
        raise RuntimeError("Single-observation deletion count gate failed")
    return keep, new_counts, index


def exact_refit(u: np.ndarray, counts: np.ndarray, initial_f: np.ndarray) -> dict:
    """Refit MBAR using all supplied native reduced biases and a warm start.

    Args:
        u: All-state reduced biases evaluated at retained observations.
        counts: Retained state counts after deletion.
        initial_f: Baseline state free energies.

    Returns:
        Normalized weights, free energies, gauge residual and overlap diagnostics.
    """
    if (
        u.shape != (len(counts), counts.sum())
        or not np.isfinite(u).all()
        or np.any(counts <= 0)
        or initial_f.shape != counts.shape
    ):
        raise ValueError("Invalid exact-bias refit inputs")
    fit = MBAR(
        u,
        counts,
        initial_f_k=initial_f,
        solver_protocol="robust",
        relative_tolerance=1e-10,
        maximum_iterations=10000,
    )
    lw = -logsumexp(np.log(counts)[:, None] + fit.f_k[:, None] - u, axis=0)
    lw -= logsumexp(lw)
    implied = -logsumexp(lw[None, :] - u, axis=1)
    residual = float(np.max(abs((implied - implied[0]) - (fit.f_k - fit.f_k[0]))))
    if residual > 1e-7:
        raise RuntimeError(f"MBAR gauge-aligned fixed-point residual {residual}")
    overlap = fit.compute_overlap()
    matrix = np.asarray(overlap["matrix"])
    components, _ = connected_components((matrix + matrix.T) > 1e-12, directed=False)
    if components != 1:
        raise RuntimeError(f"Disconnected refit overlap graph: {components}")
    eigenvalues = np.asarray(overlap["eigenvalues"])
    imaginary_error = float(np.max(abs(eigenvalues.imag)))
    if imaginary_error > 1e-10:
        raise RuntimeError(f"Non-real overlap spectrum: {imaginary_error}")
    eigenvalues = eigenvalues.real
    return dict(
        logweights=lw,
        f_k=fit.f_k,
        fixed_point_residual_kBT=residual,
        overlap_matrix=matrix,
        connected_components_at_1e_minus12=int(components),
        overlap_scalar=float(1 - eigenvalues[1]),
        overlap_eigenvalue_max_imaginary=imaginary_error,
        overlap_eigenvalues=eigenvalues.tolist(),
    )


def profile(
    s: np.ndarray,
    lw: np.ndarray,
    edges: np.ndarray,
    reference: int,
    kt: float,
    original_indices: np.ndarray,
) -> dict:
    """Compute fixed-reference PMF and per-bin support and weight concentration.

    Args:
        s: Original-coordinate projections for the selected observations.
        lw: Matching logarithmic observation weights.
        edges: Fixed baseline histogram edges.
        reference: Supported baseline reference-bin index.
        kt: Thermal energy in kcal/mol.
        original_indices: Original pooled observation IDs retained in this fit.

    Returns:
        PMF, normalized bin masses, counts, Kish ESS and dominant original IDs.
    """
    if (
        lw.shape != s.shape
        or original_indices.shape != s.shape
        or not np.isfinite(s).all()
        or not np.isfinite(lw).all()
        or not np.isfinite(edges).all()
        or np.any(np.diff(edges) <= 0)
        or s.min() < edges[0]
        or s.max() >= edges[-1]
        or not 0 <= reference < len(edges) - 1
        or not np.isfinite(kt)
        or kt <= 0
    ):
        raise ValueError("Invalid fixed-reference histogram inputs")
    lw = lw - logsumexp(lw)
    ids = np.searchsorted(edges, s, side="right") - 1
    size = len(edges) - 1
    lm = np.full(size, -np.inf)
    ess, fraction = np.zeros(size), np.zeros(size)
    dominant = np.full(size, -1, dtype=int)
    for i in range(size):
        local = np.flatnonzero(ids == i)
        if len(local):
            lm[i] = logsumexp(lw[local])
            ess[i] = np.exp(2 * lm[i] - logsumexp(2 * lw[local]))
            top = local[np.argmax(lw[local])]
            fraction[i] = np.exp(lw[top] - lm[i])
            dominant[i] = original_indices[top]
    if not np.isfinite(lm[reference]):
        raise ValueError("Frozen reference bin lost support")
    density = lm - np.log(np.diff(edges))
    pmf = -kt * (density - density[reference])
    pmf[~np.isfinite(density)] = np.nan
    return dict(
        pmf=pmf,
        mass=np.exp(lm),
        counts=np.bincount(ids, minlength=size),
        ess=ess,
        max_fraction=fraction,
        dominant_original_index=dominant,
    )
