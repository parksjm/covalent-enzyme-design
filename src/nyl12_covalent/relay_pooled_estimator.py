"""Paired projections from one global MBAR fit per independent-lineage draw."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import importlib.util
from multiprocessing import get_context

import numpy as np
from pymbar import MBAR
from scipy.special import logsumexp

from .neb_umbrella_analysis import components
from .relay_pooled_sources import FIXED
from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from .relay_umbrella_projection import project_path

KT = 0.00198720425864083 * 300
_HELPER = FIXED / "projection_sensitivity_20261126/projection.py"
_spec = importlib.util.spec_from_file_location("validated_relay_projection", _HELPER)
if _spec is None or _spec.loader is None:
    raise ImportError(_HELPER)
projection = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(projection)
_WORK: tuple | None = None


def nullable(x: np.ndarray) -> list:
    """Serialize unsupported cells explicitly, never as zero uncertainty."""
    return np.where(np.isfinite(x), x, None).tolist()


def make_projections(q: np.ndarray, config: dict) -> dict:
    """Use frozenphysicalgrids/BGPparameters;initialbranchconditioningisexplicit."""
    path = np.array(config["path_3cv_A"])
    nodes, _ = projection.uniform_nodes(path, 153)
    bgp, _ = projection.smooth_progress(q, nodes, 56.36835607868901)
    hard = project_path(q, path)["s"]
    all_frames = np.ones(len(q), bool)
    phase = (q[:, 0] >= -2) & (q[:, 0] < -0.8)
    ex = np.arange(-2.5, 3.00001, 0.1)
    ey = np.linspace(-2, -0.8, 13)
    ix = np.searchsorted(ex, q[:, 1], side="right") - 1
    iy = np.searchsorted(ey, q[:, 0], side="right") - 1
    ny = len(ey) - 1
    flat = ix * ny + iy + 0.5
    refx = int(np.searchsorted(ex, path[0, 1], side="right") - 1)
    refy = int(np.searchsorted(ey, path[0, 0], side="right") - 1)
    result = dict(
        hard=dict(
            values=hard,
            edges=np.array(config["s_bin_edges"]),
            mask=all_frames,
            reference=int(
                np.argmin(
                    abs(
                        (
                            np.array(config["s_bin_edges"])[:-1]
                            + np.array(config["s_bin_edges"])[1:]
                        )
                        / 2
                        - 0.0125
                    )
                )
            ),
        ),
        BGP=dict(values=bgp, edges=np.linspace(0, 1, 41), mask=all_frames, reference=0),
        physical_q_act=dict(values=q[:, 1], edges=ex, mask=all_frames, reference=refx),
        initial_q_act=dict(values=q[:, 1], edges=ex, mask=phase, reference=refx),
        initial_2D=dict(
            values=flat,
            edges=np.arange((len(ex) - 1) * ny + 1),
            mask=phase,
            reference=refx * ny + refy,
            x_edges=ex,
            y_edges=ey,
            shape=(len(ex) - 1, ny),
        ),
    )
    for name, p in result.items():
        values = p["values"][p["mask"]]
        if (
            len(values) == 0
            or values.min() < p["edges"][0]
            or values.max() >= p["edges"][-1]
        ):
            raise ValueError("Unsupportedorfrozen-grid-overflowprojection:" + name)
    return result


def fit_logweights(reduced: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, dict]:
    """Fit all biased states jointly; equal biases retain independent counts."""
    fit = MBAR(reduced, counts, verbose=False, solver_protocol="robust")
    overlap = fit.compute_overlap()["matrix"]
    lw = -logsumexp(np.log(counts)[:, None] + fit.f_k[:, None] - reduced, axis=0)
    lw -= logsumexp(lw)
    return lw, dict(
        f_k=fit.f_k, overlap=overlap, connected=len(components(overlap)) == 1
    )


def evaluate(
    lw: np.ndarray, projections: dict, indices: np.ndarray, support_stats: bool = True
) -> dict:
    """Rebin one fitted/resampled global weight set identically across projections."""
    result = {}
    for name, p in projections.items():
        use = p["mask"][indices]
        h = projection.local_histogram(p["values"][indices][use], lw[use], p["edges"])
        ld = np.array(h["logdensity"], float)
        ref = ld[p["reference"]]
        values = -KT * (ld - ref) if np.isfinite(ref) else np.full(len(ld), np.nan)
        result[name] = dict(
            pmf=values, ESS=h["ess"], counts=h["counts"], mass=h["mass"]
        )
        if support_stats:
            bins = (
                np.searchsorted(p["edges"], p["values"][indices][use], side="right") - 1
            )
            selected_lw = lw[use]
            result[name]["maximum_frame_weight_fraction"] = [
                None
                if not n
                else float(
                    np.exp(
                        selected_lw[bins == i].max() - logsumexp(selected_lw[bins == i])
                    )
                )
                for i, n in enumerate(h["counts"])
            ]
    return result


def initialize(reduced: np.ndarray, counts: np.ndarray, projections: dict) -> None:
    """Load large immutable arrays once per bootstrapworker, not every task."""
    global _WORK
    _WORK = (reduced, counts, projections)


def resampled_indices(counts: np.ndarray, chunk: int, seed: int) -> np.ndarray:
    """Resample whole contiguous chunks separately in every independent trajectory."""
    if np.any(counts % chunk) or np.any(counts < 2 * chunk):
        raise ValueError("RequireatleasttwoexactchunksEVERYindependentlineage")
    rng = np.random.default_rng(seed)
    chosen = []
    offset = 0
    for count in counts:
        chunks = int(count) // chunk
        picked = rng.integers(chunks, size=chunks)
        chosen.append((offset + picked[:, None] * chunk + np.arange(chunk)).ravel())
        offset += int(count)
    return np.concatenate(chosen)


def draw(task: tuple[int, int]) -> tuple[dict, bool]:
    """Return all paired projected PMFs from exactly one resampled global fit."""
    if _WORK is None:
        raise RuntimeError("Bootstrapworkernotinitialized")
    reduced, counts, projections = _WORK
    chunk, seed = task
    ids = resampled_indices(counts, chunk, seed)
    lw, fit = fit_logweights(reduced[:, ids], counts)
    estimates = evaluate(lw, projections, ids, support_stats=False)
    return {
        k: v["pmf"] if fit["connected"] else np.full_like(v["pmf"], np.nan)
        for k, v in estimates.items()
    }, fit["connected"]


def paired_bootstrap(
    reduced: np.ndarray,
    counts: np.ndarray,
    projections: dict,
    chunk: int,
    replicates: int = 200,
    workers: int = 4,
) -> dict:
    """Use existing contiguous-block protocol, with paired physical/BGP outputs.

    Args:
        reduced: K-by-N reduced native umbrella cross-energies.
        counts: Separate retained sample counts for each original seed lineage.
        projections: Frozen coordinate/bin/reference definitions.
        chunk:100frames(5ps)or50frames(2.5ps).
        replicates: Exactly200full globalrefits.
        workers: Between1and4CPUworkers.

    Returns:
        Pointwise conditional95% intervals only where every replicate has support.
    """
    if replicates != 200 or not 1 <= workers <= 4 or chunk not in (50, 100):
        raise ValueError("Unapprovedbootstrapprotocol")
    resampled_indices(counts, chunk, 0)
    tasks = [
        (chunk, int(s))
        for s in np.random.SeedSequence(20261127).generate_state(replicates)
    ]
    if workers == 1:
        initialize(reduced, counts, projections)
        draws = [draw(t) for t in tasks]
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
        values = np.asarray([d[0][name] for d in draws])
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
            scope="Pointwiseconditional95%;notsimultaneousorconvergenceproof",
        )
    return result


def point(data: dict) -> tuple[dict, dict, np.ndarray, np.ndarray]:
    """Perform the sole pooled point fit and retain offsets/weights for reuse."""
    q = data["q"]
    counts = data["counts"]
    projections = make_projections(q, data["config"])
    reduced = (
        cross_bias_energies(
            dict(zip(("q_acyl", "q_act", "q_del"), q.T, strict=True)),
            [RelayBiasWindow(**s["bias"]) for s in data["states"]],
        )
        / KT
    )
    lw, fit = fit_logweights(reduced, counts)
    result = evaluate(lw, projections, np.arange(len(q)))
    for v in result.values():
        v["pmf"] = nullable(v["pmf"])
    return (
        dict(
            projections=result,
            connected=fit["connected"],
            overlap=fit["overlap"].tolist(),
            f_k=fit["f_k"].tolist(),
            states=len(counts),
            frames=len(q),
            counts=counts.tolist(),
            CI_gate="POINTONLY-one5pschunk/newlineage"
            if not data["bootstrap_eligible"]
            else "eligible-twoormore5pschunks/lineage",
            primary_replacement=False,
            sampling_converged=False,
        ),
        projections,
        reduced,
        lw,
    )
