"""Finite-bin WHAM equations solved by the existing PyMBAR numerical solver."""

from __future__ import annotations

import numpy as np
from pymbar import MBAR
from scipy.special import logsumexp


def quantize(q: np.ndarray, width: float, shift: float = 0.0) -> tuple:
    """Return occupied cell centers, frame assignments and occupancies.

    Args:
        q: Finite physical CV values, shape (frames, dimensions).
        width: Positive cell width in Angstrom.
        shift: Grid-edge displacement in fractions of the cell width.

    Returns:
        Occupied centers, inverse frame-to-cell map and cell counts.
    """
    q = np.asarray(q, float)
    if q.ndim != 2 or not len(q) or not np.isfinite(q).all():
        raise ValueError("Finite nonempty 2D coordinates required")
    if not np.isfinite(width) or width <= 0 or not np.isfinite(shift):
        raise ValueError("Positive finite width and finite shift required")
    index = np.floor(q / width - shift).astype(np.int64)
    unique, inverse, counts = np.unique(
        index, axis=0, return_inverse=True, return_counts=True
    )
    return (unique + shift + 0.5) * width, inverse, counts


def solve_wham(
    cell_u: np.ndarray,
    inverse: np.ndarray,
    counts: np.ndarray,
    initial_f: np.ndarray | None = None,
) -> dict:
    """Solve finite-bin WHAM by expanding occupied cell energies per frame.

    Args:
        cell_u: State-by-occupied-cell reduced native biases.
        inverse: Cell index for each original frame in original state order.
        counts: Original retained frame counts for each independent history.
        initial_f: Optional warm-start state free energies in reduced units.

    Returns:
        Normalized per-frame log weights, cell masses, free energies and the
        gauge-aligned WHAM fixed-point residual in reduced units.
    """
    cell_u = np.asarray(cell_u, float)
    inverse = np.asarray(inverse, int)
    counts = np.asarray(counts, int)
    if (
        cell_u.ndim != 2
        or not np.isfinite(cell_u).all()
        or counts.shape != (len(cell_u),)
        or np.any(counts <= 0)
        or counts.sum() != len(inverse)
        or inverse.ndim != 1
        or np.any(inverse < 0)
        or np.any(inverse >= cell_u.shape[1])
    ):
        raise ValueError("Invalid state counts, finite cell biases or assignments")
    occupancy = np.bincount(inverse, minlength=cell_u.shape[1])
    if np.any(occupancy == 0):
        raise ValueError("Only occupied cells may be passed")
    fit = MBAR(
        cell_u[:, inverse],
        counts,
        initial_f_k=initial_f,
        solver_protocol="robust",
        relative_tolerance=1e-10,
        maximum_iterations=10000,
    )
    denominator = logsumexp(np.log(counts)[:, None] + fit.f_k[:, None] - cell_u, axis=0)
    cell_logmass = np.log(occupancy) - denominator
    normalizer = logsumexp(cell_logmass)
    cell_logmass -= normalizer
    logweights = -denominator[inverse] - normalizer
    implied_f = -logsumexp(cell_logmass[None, :] - cell_u, axis=1)
    residual = float(np.max(abs((implied_f - implied_f[0]) - (fit.f_k - fit.f_k[0]))))
    if residual > 1e-7:
        raise RuntimeError(f"WHAM fixed-point residual {residual} exceeds 1e-7")
    return dict(
        logweights=logweights,
        cell_logmass=cell_logmass,
        f_k=fit.f_k,
        fixed_point_residual_kBT=residual,
    )


def projected_pmf(
    values: np.ndarray,
    logweights: np.ndarray,
    edges: np.ndarray,
    reference: int,
    kt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Histogram probability masses and form a PMF at one fixed reference bin.

    Args:
        values: Original-frame or occupied-center observable values.
        logweights: Matching normalized log frame weights or log cell masses.
        edges: Frozen observable bin edges, covering every value.
        reference: Frozen supported reference-bin index.
        kt: Thermal energy in kcal/mol.

    Returns:
        Relative PMF in kcal/mol, and normalized probability masses.
    """
    values, logweights, edges = map(np.asarray, (values, logweights, edges))
    if (
        values.ndim != 1
        or logweights.shape != values.shape
        or not np.isfinite(values).all()
        or not np.isfinite(logweights).all()
        or edges.ndim != 1
        or len(edges) < 2
        or np.any(np.diff(edges) <= 0)
        or not np.isfinite(edges).all()
        or not np.isfinite(kt)
        or kt <= 0
        or not 0 <= reference < len(edges) - 1
        or values.min() < edges[0]
        or values.max() >= edges[-1]
    ):
        raise ValueError("Invalid histogram values, weights, frozen edges or reference")
    ids = np.searchsorted(edges, values, side="right") - 1
    logmass = np.array([logsumexp(logweights[ids == i]) for i in range(len(edges) - 1)])
    density = logmass - np.log(np.diff(edges))
    if not np.isfinite(density[reference]):
        raise ValueError("Reference bin has no support")
    pmf = -kt * (density - density[reference])
    pmf[~np.isfinite(density)] = np.nan
    return pmf, np.exp(logmass)
