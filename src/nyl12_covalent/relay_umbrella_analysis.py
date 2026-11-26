"""Joint two-chart MBAR and uncensored path-projected probability densities."""

from __future__ import annotations

from typing import Any

import numpy as np
from pymbar import MBAR, timeseries
from scipy.special import logsumexp

from .neb_umbrella_analysis import components
from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from .relay_umbrella_projection import path_histogram, project_path

KB = 0.00198720425864083
FEATURES = ("q_acyl_A", "q_act_A", "q_del_A")


def select_uncorrelated(values: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    """Discard detected equilibration and conservatively thin every observable.

    Args:
        values: Finite feature and orthogonal-CN time series, columns by feature.

    Returns:
        Retained indices and explicit discard/inefficiency diagnostics.
    """
    values = np.asarray(values, float)
    if values.ndim != 2 or len(values) < 4 or not np.isfinite(values).all():
        raise ValueError("At least four finite frames per accepted block required")
    estimates = [
        timeseries.detect_equilibration(column)[:2]
        for column in values.T
        if np.ptp(column) > 1e-12
    ]
    discard = max((int(item[0]) for item in estimates), default=0)
    post = values[discard:]
    inefficiency = max(
        (
            float(timeseries.statistical_inefficiency(column))
            for column in post.T
            if np.ptp(column) > 1e-12
        ),
        default=1.0,
    )
    indices = discard + np.asarray(
        timeseries.subsample_correlated_data(
            post[:, 0], g=inefficiency, conservative=True
        ),
        int,
    )
    if len(indices) < 2:
        raise ValueError("Fewer than two decorrelated post-equilibration samples")
    return indices, dict(
        raw_frames=len(values),
        discarded_frames=discard,
        statistical_inefficiency=inefficiency,
        selected_frames=len(indices),
    )


def analyze(
    windows: list[dict[str, Any]],
    blocks: list[dict[str, Any]],
    *,
    path: np.ndarray,
    temperature: float,
    edges: np.ndarray,
    es_bounds: tuple[float, float],
    acyl_bounds: tuple[float, float],
    neighbor_pairs: list[list[int]],
    minimum_duration_ps: float = 5.0,
) -> dict[str, Any]:
    """Analyze accepted compact blocks jointly, without chemical-state censoring.

    Args:
        windows: Ordered labels and ``bias`` RelayBiasWindow keyword mappings.
        blocks: Accepted records with window, values (three CVs then CN/chemical
            diagnostics), native_bias_kcal_mol, and sampled_duration_ps.
        path: Accepted three-CV polyline, in Angstrom.
        temperature: Exact common physical bath temperature.
        edges: Fixed explicit raw-s histogram edges.
        es_bounds: Declared E-S reference search interval.
        acyl_bounds: Declared acylenzyme reference search interval.
        neighbor_pairs: Explicit neighboring window indices for overlap gates.
        minimum_duration_ps: Initial full-coverage target per window.

    Returns:
        JSON-serializable coverage, eight gates, overlap and preliminary PMF.
    """
    if not np.isfinite(temperature) or temperature <= 0 or len(windows) < 2:
        raise ValueError("Positive temperature and at least two windows required")
    labels = [w["label"] for w in windows]
    if len(set(labels)) != len(labels):
        raise ValueError("Window labels must be unique")
    project_path(np.empty((0, 3)), path)
    edges = np.asarray(edges, float)
    if (
        edges.ndim != 1
        or len(edges) < 3
        or not np.isfinite(edges).all()
        or np.any(np.diff(edges) <= 0)
    ):
        raise ValueError("Declare increasing finite path histogram edges")
    if not es_bounds[0] < es_bounds[1] or not acyl_bounds[0] < acyl_bounds[1]:
        raise ValueError("Declare nonempty E-S and acylenzyme search intervals")
    if any(
        len(pair) != 2
        or pair[0] == pair[1]
        or any(type(index) is not int or not 0 <= index < len(labels) for index in pair)
        for pair in neighbor_pairs
    ):
        raise ValueError("Neighbor pairs must identify two distinct manifest states")
    biases = [RelayBiasWindow(**w["bias"]) for w in windows]
    grouped: dict[str, list[np.ndarray]] = {label: [] for label in labels}
    durations = dict.fromkeys(labels, 0.0)
    summaries = []
    for block in blocks:
        label = block["window"]
        if label not in grouped:
            raise ValueError("Accepted block names an unknown window")
        values = np.asarray(block["values"], float)
        native = np.asarray(block["native_bias_kcal_mol"], float)
        energy = cross_bias_energies(
            dict(zip(("q_acyl", "q_act", "q_del"), values[:, :3].T, strict=True)),
            [biases[labels.index(label)]],
        )[0]
        mask = np.asarray(
            block.get("native_bias_evaluated", np.ones(len(values), dtype=bool)), bool
        )
        if mask.shape != energy.shape or not mask[:-1].all() or not mask.any():
            raise ValueError(
                "Native parity can be unevaluated only for the terminal frame"
            )
        eps = float(block.get("native_cv_quantization_A", 0.0)) + float(
            block.get("coordinate_serialization_tolerance_A", 0.0)
        )
        if not np.isfinite(eps) or eps < 0:
            raise ValueError(
                "Native CV quantization/serialization bound must be finite and nonnegative"
            )
        allowed = np.zeros_like(energy)
        own_bias = biases[labels.index(label)]
        for name, center, force in zip(
            own_bias.feature_names,
            own_bias.centers_A,
            own_bias.force_constants,
            strict=True,
        ):
            column = values[:, ("q_acyl", "q_act", "q_del").index(name)]
            allowed += (
                2
                * force
                * np.minimum(np.abs(column - center) + eps, own_bias.outer_extent_A)
                * eps
            )
        energy_rounding = float(block.get("native_energy_quantization_kcal_mol", 0.0))
        if not np.isfinite(energy_rounding) or energy_rounding < 0:
            raise ValueError(
                "Native energy rounding bound must be finite and nonnegative"
            )
        allowed += energy_rounding
        numerical = 32 * np.finfo(float).eps * np.maximum(1, np.abs(energy))
        if (
            native.shape != energy.shape
            or not np.isfinite(native[mask]).all()
            or np.any(
                np.abs(native[mask] - energy[mask]) > allowed[mask] + numerical[mask]
            )
        ):
            raise ValueError("Native Amber/computed per-frame bias parity failed")
        if "native_cv_values_A" in block:
            native_cv = np.asarray(block["native_cv_values_A"], float)
            actual = np.column_stack(
                [
                    values[:, ("q_acyl", "q_act", "q_del").index(name)]
                    for name in own_bias.feature_names
                ]
            )
            if (
                native_cv.shape != actual.shape
                or not np.isfinite(native_cv[mask]).all()
                or np.any(np.abs(native_cv[mask] - actual[mask]) > eps + 1e-12)
            ):
                raise ValueError("Independent native CV parity failed")
        duration = float(block["sampled_duration_ps"])
        if not np.isfinite(duration) or duration <= 0:
            raise ValueError("Positive finite accepted duration required")
        indices, summary = select_uncorrelated(values)
        grouped[label].append(values[indices])
        durations[label] += duration
        summaries.append(
            dict(
                window=label,
                native_bias_evaluated_frames=int(mask.sum()),
                native_bias_unevaluated_terminal_frames=int((~mask).sum()),
                native_bias_source=block.get(
                    "native_bias_source", "exact_native_energy"
                ),
                native_bias_max_error_kcal_mol=float(
                    np.max(np.abs(native[mask] - energy[mask]))
                ),
                native_bias_max_bound_kcal_mol=float(np.max(allowed[mask])),
                **summary,
            )
        )
    counts = np.array([sum(len(part) for part in grouped[label]) for label in labels])
    covered = all(durations[label] >= minimum_duration_ps for label in labels)
    gates = {
        name: dict(passed=None, reason="not evaluated")
        for name in (
            "neighbor_histogram_overlap",
            "connected_mbar",
            "effective_samples",
            "block_stability",
            "bootstrap_95CI",
            "WHAM_agreement",
            "independent_starts",
            "orthogonal_stationarity",
        )
    }
    gates["effective_samples"] = dict(
        passed=bool(np.all(counts >= 100)), minimum=100, N_k=counts.tolist()
    )
    result = dict(
        status="readiness_no_accepted_data" if not blocks else "incomplete_coverage",
        full_5ps_coverage=covered,
        duration_ps=durations,
        N_k=counts.tolist(),
        blocks=summaries,
        gates=gates,
        pmf=None,
        converged=False,
        scope="Conditional modeled bound reaction channel; not assay kcat",
        projection="Equal-Angstrom nearest 3CV polyline; terminal tangents extended",
        path_length_A=float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum()),
        retained_chemical_branches="All supplied frames; no proton-owner censoring",
    )
    if not blocks or not covered:
        return result
    samples = np.concatenate([np.concatenate(grouped[label]) for label in labels])
    raw = dict(zip(("q_acyl", "q_act", "q_del"), samples[:, :3].T, strict=True))
    reduced = cross_bias_energies(raw, biases) / (KB * temperature)
    mbar = MBAR(reduced, counts, verbose=False, solver_protocol="robust")
    overlap = mbar.compute_overlap()["matrix"]
    groups = components(overlap)
    gates["connected_mbar"] = dict(passed=len(groups) == 1, edge_threshold=0.001)
    result.update(
        overlap=overlap.tolist(), components=groups, status="disconnected_no_global_PMF"
    )
    # Empirical histogram intersection, with one common 3CV grid for all charts.
    grid = [
        np.linspace(samples[:, j].min() - 1e-8, samples[:, j].max() + 1e-8, 13)
        for j in range(3)
    ]
    partitions = np.split(samples[:, :3], np.cumsum(counts)[:-1])
    hist = [np.histogramdd(part, bins=grid)[0] / len(part) for part in partitions]
    intersections = [
        float(np.minimum(hist[i], hist[j]).sum()) for i, j in neighbor_pairs
    ]
    gates["neighbor_histogram_overlap"] = dict(
        passed=bool(intersections and min(intersections) >= 0.2),
        minimum=0.2,
        pairs=neighbor_pairs,
        intersections=intersections,
        grid_edges=[g.tolist() for g in grid],
    )
    if len(groups) != 1:
        return result
    logw = -logsumexp(np.log(counts)[:, None] + mbar.f_k[:, None] - reduced, axis=0)
    weights = np.exp(logw - logsumexp(logw))
    projection = project_path(samples[:, :3], path)
    result.update(
        status="preliminary_connected_not_converged",
        unbiased_ESS=float(1 / np.square(weights).sum()),
        pmf=path_histogram(
            projection["s"], weights, edges, temperature, es_bounds, acyl_bounds
        ),
        projection_diagnostics={
            key: [float(x) if np.isfinite(x) else None for x in value]
            for key, value in projection.items()
        },
        diagnostic_columns=samples[:, 3:].tolist(),
        common_U0="Cancelled only for exactly matched physical Hamiltonian and temperature",
        bias_free_energies_are_not_PMF=True,
    )
    return result
