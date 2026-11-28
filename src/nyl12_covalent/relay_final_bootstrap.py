"""Isolated count-preserving circular moving blocks for paired MBAR projections."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np

from .relay_pooled_estimator import evaluate, fit_logweights, nullable, paired_bootstrap

_WORK: tuple | None = None


def circular_indices(counts: np.ndarray, block_frames: int, seed: int) -> np.ndarray:
    """Draw circular contiguous blocks independently per lineage, retaining exact counts.

    Args:
        counts: Actual independent trajectory sample counts, not merged bias states.
        block_frames:100for5psor200for10psat0.05ps/sample.
        seed: Explicit reproducible draw seed.

    Returns:
        Exactlysum(counts)indices. No original tail is removed; all samples canstartblocks.
    """
    counts = np.asarray(counts, int)
    if (
        counts.ndim != 1
        or block_frames not in (100, 200)
        or np.any(counts < 2 * block_frames)
    ):
        raise ValueError("Require>=twocircularblocksperindependentlineage")
    rng = np.random.default_rng(seed)
    result = []
    offset = 0
    for count in counts:
        starts = rng.integers(0, int(count), size=int(np.ceil(count / block_frames)))
        local = ((starts[:, None] + np.arange(block_frames)) % int(count)).ravel()[
            : int(count)
        ]
        result.append(local + offset)
        offset += int(count)
    return np.concatenate(result)


def initialize(reduced: np.ndarray, counts: np.ndarray, projections: dict) -> None:
    """Load large immutable fit arrays only once per worker."""
    global _WORK
    _WORK = (reduced, counts, projections)


def draw(task: tuple[int, int]) -> tuple[dict, bool]:
    """Fit one count-preserving draw and use its global weights for all coordinates."""
    if _WORK is None:
        raise RuntimeError("Finalbootstrapworkernotinitialized")
    reduced, counts, projections = _WORK
    block, seed = task
    ids = circular_indices(counts, block, seed)
    lw, fit = fit_logweights(reduced[:, ids], counts)
    estimates = evaluate(lw, projections, ids, support_stats=False)
    return {
        name: value["pmf"] if fit["connected"] else np.full_like(value["pmf"], np.nan)
        for name, value in estimates.items()
    }, fit["connected"]


def circular_bootstrap(
    reduced: np.ndarray,
    counts: np.ndarray,
    projections: dict,
    block_frames: int,
    workers: int = 4,
    replicates: int = 200,
) -> dict:
    """Compute paired conditional intervals with arbitrary retained trajectory lengths.

    Args:
        reduced: Native umbrella reduced crossenergies.
        counts: Independent lineage counts, including1100/500frames.
        projections: Same frozen coordinate/grids for every resample.
        block_frames:100or200samples, never rounded down todividecounts.
        workers: AtmostfourCPUworkers.
        replicates: Exactly200fullglobalMBARrefits.

    Returns:
        Paired pointwiseconditional95% intervals; circularassumption anddrift remainexplicit.
    """
    if not 1 <= workers <= 4 or replicates != 200:
        raise ValueError("Require200refitsandmax4CPU")
    circular_indices(counts, block_frames, 0)
    tasks = [
        (block_frames, int(s))
        for s in np.random.SeedSequence(20261127).generate_state(replicates)
    ]
    if workers == 1:
        initialize(reduced, counts, projections)
        draws = [draw(task) for task in tasks]
    else:
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=get_context("spawn"),
            initializer=initialize,
            initargs=(reduced, counts, projections),
        ) as pool:
            draws = list(pool.map(draw, tasks))
    result = {}
    for name in projections:
        values = np.array([d[0][name] for d in draws])
        support = np.isfinite(values).mean(axis=0)
        lower = np.full(values.shape[1], np.nan)
        upper = lower.copy()
        full = support == 1
        if full.any():
            lower[full], upper[full] = np.percentile(
                values[:, full], [2.5, 97.5], axis=0
            )
        result[name] = dict(
            lower_95=nullable(lower),
            upper_95=nullable(upper),
            support_fraction=support.tolist(),
            replicates_pmf=nullable(values),
            disconnected_replicates=sum(not d[1] for d in draws),
        )
    return dict(
        method="Circularmovingblockbootstrap;independentlineages;countpreserving",
        block_frames=block_frames,
        block_ps=block_frames * 0.05,
        replicates=replicates,
        workers=workers,
        counts=counts.tolist(),
        last5ps_dropped=False,
        projections=result,
        interpretation="Pointwiseconditional95%;circularboundaryassumptiondoesnotrepairdriftorproveconvergence",
    )


def nonoverlap_baseline(
    reduced: np.ndarray, counts: np.ndarray, projections: dict, workers: int = 4
) -> dict:
    """Keep the existing5psnonoverlapping protocol as an explicitly distinct baseline."""
    return dict(
        method="Existingindependentlineagenonoverlapping5psblockbaseline",
        projections=paired_bootstrap(
            reduced, counts, projections, 100, workers=workers
        ),
        interpretation="Comparatortocircular5/10ps;NOTindependent-samplingorstationarityproof",
    )
