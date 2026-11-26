"""Refit MBAR after resampling completed five-ps blocks within each window."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from typing import Any

import numpy as np
from pymbar import MBAR
from scipy.special import logsumexp

from .relay_umbrella_analysis import KB


def _refit(arguments: tuple[Any, ...]) -> list[float]:
    reduced, s, edges, temperature, reference_bin, groups, seed = arguments
    rng = np.random.default_rng(seed)
    indices, counts = [], []
    for blocks in groups:
        selected = np.concatenate(
            [blocks[i] for i in rng.integers(len(blocks), size=len(blocks))]
        )
        indices.append(selected)
        counts.append(len(selected))
    chosen = np.concatenate(indices)
    potentials = reduced[:, chosen]
    mbar = MBAR(potentials, np.array(counts), verbose=False, solver_protocol="robust")
    logweights = -logsumexp(
        np.log(counts)[:, None] + mbar.f_k[:, None] - potentials, axis=0
    )
    bins = np.searchsorted(edges, s[chosen], side="right") - 1
    values = np.full(len(edges) - 1, np.nan)
    for i in np.unique(bins):
        values[i] = (
            -KB
            * temperature
            * (logsumexp(logweights[bins == i]) - np.log(np.diff(edges)[i]))
        )
    return (values - values[reference_bin]).tolist()


def block_bootstrap(
    reduced: np.ndarray,
    counts: np.ndarray,
    s: np.ndarray,
    edges: np.ndarray,
    temperature: float,
    reference_s: float,
    block_lengths: list[list[int]],
    *,
    replicates: int,
    seed: int,
    workers: int,
    resampling_unit: str = "completed_5ps_block",
) -> dict[str, Any]:
    """Calculate pointwise percentile intervals using complete block refits.

    Args:
        reduced: Cross-evaluated reduced biases, K by N, in window order.
        counts: Original per-window selected sample counts.
        s: Original selected path projections.
        edges: Fixed histogram edges for every replicate.
        temperature: Physical bath temperature in kelvin.
        reference_s: Fixed reference-bin center, never re-minimized.
        block_lengths: Selected frame counts per completed block in each window.
        replicates: Number of independent bootstrap refits, at least 50.
        seed: Reproducible NumPy seed.
        workers: Number of CPU processes, explicitly bounded by caller.
        resampling_unit: Completed five-ps blocks or individually decorrelated
            configurations; the latter requires singleton groups and assumes
            the upstream whole-window correlation estimate is adequate.

    Returns:
        Pointwise conditional percentile intervals, per-bin support fractions,
        and full replicates. Missing-support draws are retained as null rather
        than silently removed from the support assessment. At short durations,
        four blocks/window provide limited evidence about slow relaxation.
    """
    if replicates < 50 or not 1 <= workers <= 8:
        raise ValueError("At least50 replicates and1–8 workers required")
    if resampling_unit not in {"completed_5ps_block", "decorrelated_configuration"}:
        raise ValueError("Unknown bootstrap resampling unit")
    if resampling_unit == "decorrelated_configuration" and any(
        length != 1 for lengths in block_lengths for length in lengths
    ):
        raise ValueError(
            "Decorrelated-configuration bootstrap requires singleton groups"
        )
    counts = np.asarray(counts, int)
    reduced, s, edges = map(lambda value: np.asarray(value, float), (reduced, s, edges))
    if (
        reduced.shape != (len(counts), len(s))
        or counts.sum() != len(s)
        or len(block_lengths) != len(counts)
        or np.any(counts <= 0)
        or not np.isfinite(reduced).all()
        or not np.isfinite(s).all()
        or not np.isfinite(temperature)
        or temperature <= 0
    ):
        raise ValueError(
            "Finite matching samples, block counts and temperature required"
        )
    if (
        edges.ndim != 1
        or len(edges) < 3
        or np.any(np.diff(edges) <= 0)
        or not np.isfinite(edges).all()
        or s.min() < edges[0]
        or s.max() >= edges[-1]
    ):
        raise ValueError("Finite increasing edges must cover all samples")
    centers = (edges[:-1] + edges[1:]) / 2
    reference = np.flatnonzero(np.isclose(centers, reference_s, rtol=0, atol=1e-12))
    if len(reference) != 1:
        raise ValueError("Reference must identify exactly one bin center")
    groups, offset = [], 0
    for count, lengths in zip(counts, block_lengths, strict=True):
        if len(lengths) < 2 or min(lengths) <= 0 or sum(lengths) != count:
            raise ValueError("At least two nonempty blocks must partition each window")
        window = []
        for length in lengths:
            window.append(np.arange(offset, offset + length))
            offset += length
        groups.append(window)
    seeds = np.random.SeedSequence(seed).generate_state(replicates)
    arguments = [
        (reduced, s, edges, temperature, reference[0], groups, int(state))
        for state in seeds
    ]
    if workers == 1:
        values = np.asarray([_refit(args) for args in arguments])
    else:
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=get_context("spawn")
        ) as pool:
            values = np.asarray(list(pool.map(_refit, arguments)))
    supported = np.isfinite(values)
    support = supported.mean(axis=0)
    lower, upper = np.full(len(centers), np.nan), np.full(len(centers), np.nan)
    for index in np.flatnonzero(support > 0):
        lower[index], upper[index] = np.percentile(
            values[supported[:, index], index], [2.5, 97.5]
        )

    def nullable(array: np.ndarray) -> list[float | None]:
        return [float(value) if np.isfinite(value) else None for value in array]

    return dict(
        method=(
            "Within-window completed5ps block bootstrap; fixed source thinning; MBAR refit"
            if resampling_unit == "completed_5ps_block"
            else "Window-stratified decorrelated-configuration bootstrap; MBAR refit"
        ),
        resampling_unit=resampling_unit,
        replicates=replicates,
        seed=seed,
        workers=workers,
        selected_block_lengths=block_lengths,
        lower_95_kcal_mol=nullable(lower),
        upper_95_kcal_mol=nullable(upper),
        support_fraction=support.tolist(),
        replicates_kcal_mol=[nullable(row) for row in values],
        reference_s=reference_s,
        interpretation=(
            "Percentiles conditional on bin/reference support; display intervals only"
            " where support_fraction==1. Not a simultaneous band or convergence proof."
            " Conditional on fixed burn-in and decorrelation; cannot exclude slower correlations."
        ),
    )
