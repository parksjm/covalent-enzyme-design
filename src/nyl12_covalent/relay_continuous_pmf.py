"""One burn-in and correlation estimate per continuous production window."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import numpy as np
from scipy.special import logsumexp

from .neb_umbrella_analysis import components
from .relay_pmf_block_bootstrap import block_bootstrap
from .relay_pmf_uncertainty import histogram_uncertainty
from .relay_umbrella_analysis import KB, select_uncorrelated
from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from .relay_umbrella_projection import path_histogram, project_path


def continuous_selections(
    windows: list[dict[str, Any]],
    blocks: list[dict[str, Any]],
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Join production segments in physical-time order before any selection.

    Args:
        windows: Ordered umbrella-state definitions.
        blocks: Hash-validated production blocks with original time-stamped rows.

    Returns:
        Window-ordered selected arrays, per-window counts, and one selection
        summary/window. No segmentwise equilibration or thinning is applied.
    """
    known = {row["label"] for row in windows}
    if any(block["window"] not in known for block in blocks):
        raise ValueError("Unknown window in continuous production blocks")
    selections, summaries = [], []
    for window in windows:
        parts = [block for block in blocks if block["window"] == window["label"]]
        if not parts:
            raise ValueError("Every continuous window needs accepted production data")
        parts.sort(key=lambda block: float(block["original_rows"][0]["time_ps"]))
        times = np.concatenate(
            [[float(row["time_ps"]) for row in part["original_rows"]] for part in parts]
        )
        values = np.concatenate([part["values"] for part in parts])
        spacing = np.diff(times)
        if (
            len(times) != len(values)
            or not np.isfinite(times).all()
            or np.any(spacing <= 0)
            or not np.allclose(spacing, np.median(spacing), rtol=0, atol=1e-4)
        ):
            raise ValueError("Production times have a gap, duplicate, or discontinuity")
        indices, summary = select_uncorrelated(values)
        selections.append(values[indices])
        dt = float(np.median(spacing))
        summaries.append(
            dict(
                window=window["label"],
                **summary,
                accepted_blocks=len(parts),
                duration_ps=float(sum(part["sampled_duration_ps"] for part in parts)),
                frame_interval_ps=dt,
                first_frame_time_ps=float(times[0]),
                last_frame_time_ps=float(times[-1]),
                burn_in_ps=float(summary["discarded_frames"] * dt),
                selected_indices=indices.tolist(),
                retained_duration_ps=float(
                    (len(times) - summary["discarded_frames"]) * dt
                ),
            )
        )
    counts = np.array([len(values) for values in selections], int)
    return np.concatenate(selections), counts, summaries


def analyze_continuous(
    config: dict[str, Any],
    blocks: list[dict[str, Any]],
    parent: dict[str, Any],
    *,
    replicates: int,
    workers: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Refit the same PMF after continuous-window equilibration/decorrelation.

    Args:
        config: Unchanged source request with Hamiltonian, biases, bins and path.
        blocks: Hash-validated source production records in any order.
        parent: Frozen blockwise report to compare, never modified.
        replicates: Explicit bootstrap refits of window-stratified selected frames.
        workers: Bounded CPU processes for those refits.

    Returns:
        New analysis and uncertainty reports. Analytical intervals are directly
        comparable to the parent's analytical intervals. The supplemental
        bootstrap assumes effective independence after whole-window thinning.
    """
    if len(blocks) != len(parent["accepted_sources"]) or any(
        block["receipt_sha256"] != source["receipt_sha256"]
        or block["source_sha256"] != source["source_sha256"]
        for block, source in zip(blocks, parent["accepted_sources"], strict=True)
    ):
        raise ValueError(
            "Continuous analysis must use exactly the accepted parent data"
        )
    if config["temperature_K"] != parent["temperature_K"]:
        raise ValueError("Source physical temperature changed")
    samples, counts, summaries = continuous_selections(config["windows"], blocks)
    raw = dict(zip(("q_acyl", "q_act", "q_del"), samples[:, :3].T, strict=True))
    biases = [RelayBiasWindow(**window["bias"]) for window in config["windows"]]
    reduced = cross_bias_energies(raw, biases) / (KB * config["temperature_K"])
    projection = project_path(samples[:, :3], np.array(config["path_3cv_A"]))
    reference = parent["checkpoint"]["plot_common_es_reference_s"]
    edges = np.array(config["s_bin_edges"])
    uncertainty = histogram_uncertainty(
        reduced,
        counts,
        projection["s"],
        edges,
        config["temperature_K"],
        reference,
    )
    overlap = np.array(uncertainty["overlap"])
    groups = components(overlap)
    if len(groups) != 1:
        raise ValueError("Continuous-window MBAR is disconnected; no global PMF")
    logweights = -logsumexp(
        np.log(counts)[:, None] + np.array(uncertainty["mbar_f_k"])[:, None] - reduced,
        axis=0,
    )
    weights = np.exp(logweights - logsumexp(logweights))
    pmf = path_histogram(
        projection["s"],
        weights,
        edges,
        config["temperature_K"],
        tuple(config["es_search_bounds"]),
        tuple(config["acyl_search_bounds"]),
    )
    # Cover all common raw frames, including newly selected tail configurations.
    old_gate = parent["gates"]["neighbor_histogram_overlap"]
    all_raw = np.concatenate([block["values"][:, :3] for block in blocks])
    grid = [
        np.linspace(all_raw[:, j].min() - 1e-8, all_raw[:, j].max() + 1e-8, 13)
        for j in range(3)
    ]
    partitions = np.split(samples[:, :3], np.cumsum(counts)[:-1])
    histograms = [np.histogramdd(part, bins=grid)[0] / len(part) for part in partitions]
    intersections = [
        float(np.minimum(histograms[i], histograms[j]).sum())
        for i, j in config["neighbor_pairs"]
    ]
    result = deepcopy(parent)
    for gate in result["gates"].values():
        gate.clear()
        gate.update(passed=None, reason="Not evaluated for continuous-window estimator")
    result["gates"]["effective_samples"] = dict(
        passed=bool(np.all(counts >= 100)),
        minimum=100,
        N_k=counts.tolist(),
    )
    result["gates"]["connected_mbar"] = dict(passed=True, edge_threshold=0.001)
    result["gates"]["neighbor_histogram_overlap"] = dict(
        passed=min(intersections) >= old_gate["minimum"],
        minimum=old_gate["minimum"],
        pairs=config["neighbor_pairs"],
        intersections=intersections,
        grid_edges=[axis.tolist() for axis in grid],
        grid_scope="Common raw-source support; not matched to archived selected-data grid",
    )
    result.update(
        status="preliminary_connected_not_converged",
        converged=False,
        N_k=counts.tolist(),
        overlap=overlap.tolist(),
        components=groups,
        unbiased_ESS=uncertainty["unbiased_ESS"],
        pmf=pmf,
        selection_policy="continuous_window",
        window_selection=summaries,
        diagnostic_columns=samples[:, 3:].tolist(),
        projection_diagnostics={
            key: value.tolist() for key, value in projection.items()
        },
        blocks=[
            dict(
                window=block["window"],
                raw_frames=len(block["values"]),
                scope="Native source acceptance only; selection reported per window",
            )
            for block in blocks
        ],
    )
    result.pop("analysis_revision", None)
    uncertainty.update(
        selection_policy="continuous_window",
        temperature_K=config["temperature_K"],
        N_k=counts.tolist(),
        window_selection=summaries,
        physical_hamiltonian_sha256=config["physical_hamiltonian_sha256"],
        path_3cv_sha256=parent["path_3cv_sha256"],
    )
    if replicates:
        uncertainty["window_bootstrap"] = block_bootstrap(
            reduced,
            counts,
            projection["s"],
            edges,
            config["temperature_K"],
            reference,
            [[1] * int(count) for count in counts],
            replicates=replicates,
            seed=20261126,
            workers=workers,
            resampling_unit="decorrelated_configuration",
        )
    return result, uncertainty
