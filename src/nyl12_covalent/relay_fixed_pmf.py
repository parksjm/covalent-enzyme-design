"""All-frame MBAR estimates with within-window contiguous-block uncertainty."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from typing import Any

import numpy as np
from pymbar import MBAR
from scipy.special import logsumexp

from .neb_umbrella_analysis import components
from .relay_umbrella_analysis import KB


def fit_samples(
    reduced: np.ndarray,
    counts: np.ndarray,
    s: np.ndarray,
    edges: np.ndarray,
    temperature: float,
    reference_s: float,
) -> dict[str, Any]:
    """Fit the point estimate without treating correlated frames as independent.

    Args:
        reduced: Cross-evaluated reduced umbrella biases, shape K by N.
        counts: Actual retained frame count in each biased trajectory.
        s: Projection of each retained configuration onto the unchanged path.
        edges: Fixed histogram edges strictly covering all projections.
        temperature: Common physical temperature in kelvin.
        reference_s: Fixed reference-bin center shared across every refit.

    Returns:
        Stable histogram estimates and overlap/support diagnostics. No IID
        analytical uncertainty is calculated. Correlation affects precision,
        while stationarity and ergodicity remain requirements for consistency.
        Disconnected fits are explicitly identified for bootstrap rejection.
    """
    reduced, s, edges = map(lambda x: np.asarray(x, float), (reduced, s, edges))
    counts = np.asarray(counts)
    if (
        counts.ndim != 1
        or not np.issubdtype(counts.dtype, np.integer)
        or np.any(counts <= 0)
        or counts.sum() != len(s)
        or reduced.shape != (len(counts), len(s))
        or not np.isfinite(reduced).all()
        or not np.isfinite(s).all()
        or temperature <= 0
        or not np.isfinite(temperature)
    ):
        raise ValueError(
            "Require finite matching arrays, integer counts and temperature"
        )
    if (
        edges.ndim != 1
        or len(edges) < 3
        or not np.isfinite(edges).all()
        or np.any(np.diff(edges) <= 0)
        or s.min() < edges[0]
        or s.max() >= edges[-1]
    ):
        raise ValueError("Increasing finite edges must strictly cover projections")
    centers = (edges[:-1] + edges[1:]) / 2
    reference = np.flatnonzero(np.isclose(centers, reference_s, rtol=0, atol=1e-12))
    if len(reference) != 1:
        raise ValueError("Reference must identify exactly one histogram-bin center")
    mbar = MBAR(reduced, counts, verbose=False, solver_protocol="robust")
    overlap = mbar.compute_overlap()["matrix"]
    connected = len(components(overlap)) == 1
    logw = -logsumexp(np.log(counts)[:, None] + mbar.f_k[:, None] - reduced, axis=0)
    bins = np.searchsorted(edges, s, side="right") - 1
    population = np.bincount(bins, minlength=len(centers))
    values, ess = np.full(len(centers), np.nan), np.zeros(len(centers))
    for i in np.flatnonzero(population):
        local = logw[bins == i]
        values[i] = -KB * temperature * (logsumexp(local) - np.log(np.diff(edges)[i]))
        ess[i] = np.exp(2 * logsumexp(local) - logsumexp(2 * local))
    values -= values[reference[0]]
    if not connected:
        values[:] = np.nan
    return dict(
        centers=centers.tolist(),
        edges=edges.tolist(),
        reference_s=reference_s,
        pmf_kcal_mol=nullable(values),
        frames_per_bin=population.tolist(),
        weight_ess_per_bin=ess.tolist(),
        weight_ess_scope="Kish weight concentration only; no temporal-correlation correction",
        global_weight_ess=float(np.exp(2 * logsumexp(logw) - logsumexp(2 * logw))),
        connected=connected,
        overlap=overlap.tolist(),
        N_k=counts.tolist(),
        sampling_converged=False,
        global_barrier_accepted=False,
    )


def nullable(values: np.ndarray) -> list[float | None]:
    """Serialize unsupported numerical estimates as explicit null values."""
    return [float(x) if np.isfinite(x) else None for x in values]


def _draw(arguments: tuple[Any, ...]) -> tuple[list[float | None], bool]:
    reduced, counts, s, edges, temperature, reference_s, chunk_frames, seed = arguments
    rng = np.random.default_rng(seed)
    chosen, offset = [], 0
    for count in counts:
        chunks = int(count) // chunk_frames
        picked = rng.integers(chunks, size=chunks)
        chosen.append(
            (offset + picked[:, None] * chunk_frames + np.arange(chunk_frames)).ravel()
        )
        offset += int(count)
    indices = np.concatenate(chosen)
    report = fit_samples(
        reduced[:, indices], counts, s[indices], edges, temperature, reference_s
    )
    return report["pmf_kcal_mol"], report["connected"]


def contiguous_bootstrap(
    reduced: np.ndarray,
    counts: np.ndarray,
    s: np.ndarray,
    edges: np.ndarray,
    temperature: float,
    reference_s: float,
    *,
    chunk_frames: int,
    replicates: int,
    workers: int,
    seed: int,
) -> dict[str, Any]:
    """Resample whole contiguous chunks independently within each window.

    Args:
        reduced: K-by-N cross-evaluated reduced umbrella biases.
        counts: Actual frame counts in trajectory order, without thinning.
        s: Matching path projections.
        edges: Fixed histogram edges.
        temperature: Physical temperature in kelvin.
        reference_s: Fixed common reference bin.
        chunk_frames: Exact equal chunk size dividing every trajectory.
        replicates: Number of full MBAR refits, at least 50.
        workers: Explicit process count, at most eight.
        seed: Reproducible random seed.

    Returns:
        Pointwise percentile intervals conditional on observed sampled states.
        Unsupported/disconnected draws remain null and enter support fractions.
        Whole-trajectory nonstationarity is not corrected by this procedure.
    """
    counts = np.asarray(counts, int)
    if (
        chunk_frames < 1
        or np.any(counts % chunk_frames)
        or np.any(counts < 2 * chunk_frames)
        or replicates < 50
        or not 1 <= workers <= 8
    ):
        raise ValueError("Require >=2 exact chunks/window, >=50 refits and1–8 workers")
    fit_samples(reduced, counts, s, edges, temperature, reference_s)
    seeds = np.random.SeedSequence(seed).generate_state(replicates)
    arguments = [
        (reduced, counts, s, edges, temperature, reference_s, chunk_frames, int(value))
        for value in seeds
    ]
    if workers == 1:
        draws = [_draw(args) for args in arguments]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=get_context("spawn")
        ) as pool:
            draws = list(pool.map(_draw, arguments))
    values = np.array([row[0] for row in draws], float)
    support = np.isfinite(values).mean(axis=0)
    lower, upper = np.full(values.shape[1], np.nan), np.full(values.shape[1], np.nan)
    for i in np.flatnonzero(support == 1):
        lower[i], upper[i] = np.percentile(values[:, i], [2.5, 97.5])
    return dict(
        method="Within-window nonoverlapping contiguous-block bootstrap; full MBAR refit",
        chunk_frames=chunk_frames,
        chunks_per_window=(counts // chunk_frames).tolist(),
        replicates=replicates,
        workers=workers,
        seed=seed,
        lower_95_kcal_mol=nullable(lower),
        upper_95_kcal_mol=nullable(upper),
        support_fraction=support.tolist(),
        disconnected_replicates=sum(not row[1] for row in draws),
        replicates_kcal_mol=[nullable(row) for row in values],
        interval_scope="Pointwise conditional95%; not simultaneous, not a barrier CI or convergence proof",
    )
