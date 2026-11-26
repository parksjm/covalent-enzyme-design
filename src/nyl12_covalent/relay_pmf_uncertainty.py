"""Statistical errors and local support for an unchanged relay histogram PMF."""

from __future__ import annotations

from typing import Any

import numpy as np
from pymbar import FES
from scipy.special import logsumexp

from .relay_umbrella_analysis import KB, select_uncorrelated
from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from .relay_umbrella_projection import project_path


def selected_samples(
    windows: list[dict[str, Any]], blocks: list[dict[str, Any]], report: dict[str, Any]
) -> tuple[np.ndarray, np.ndarray]:
    """Recover exactly the source estimator's selections, checking every block.

    Args:
        windows: Ordered source umbrella definitions.
        blocks: Hash-validated source blocks from ``load_block``.
        report: Frozen matching MBAR report with per-block selection summaries.

    Returns:
        Window-ordered samples and integer counts per window.
    """
    if len(blocks) != len(report["blocks"]):
        raise ValueError("Source blocks and selection summaries differ")
    grouped: dict[str, list[np.ndarray]] = {w["label"]: [] for w in windows}
    for block, expected in zip(blocks, report["blocks"], strict=True):
        indices, summary = select_uncorrelated(block["values"])
        if block["window"] != expected["window"] or any(
            summary[key] != expected[key] for key in summary
        ):
            raise ValueError("Frozen per-block decorrelation selection changed")
        grouped[block["window"]].append(block["values"][indices])
    counts = np.array([sum(map(len, parts)) for parts in grouped.values()], int)
    if np.any(counts == 0) or counts.tolist() != report["N_k"]:
        raise ValueError("Frozen window sample counts changed or are empty")
    return np.concatenate([np.concatenate(parts) for parts in grouped.values()]), counts


def histogram_uncertainty(
    reduced: np.ndarray,
    counts: np.ndarray,
    s: np.ndarray,
    edges: np.ndarray,
    temperature: float,
    reference_s: float,
) -> dict[str, Any]:
    """Estimate histogram covariance and local reweighting support with PyMBAR.

    Args:
        reduced: Cross-evaluated umbrella biases divided by kBT, shape K by N.
        counts: Decorrelated configuration counts for each sampled state.
        s: Matching N path projections, in window order.
        edges: Fixed histogram edges covering all samples, including endpoints.
        temperature: Common physical temperature in kelvin.
        reference_s: Fixed sampled reference-bin coordinate.

    Returns:
        Pointwise asymptotic standard errors and approximate normal 95% limits,
        counts, bin-local Kish weight ESS, and contributing window counts.
        Empty bins remain null. These errors are conditional on sampled states
        and the input decorrelation, not bounds on slow-mode or model errors.
    """
    reduced, counts = np.asarray(reduced, float), np.asarray(counts, int)
    s, edges = np.asarray(s, float), np.asarray(edges, float)
    if (
        reduced.shape != (len(counts), len(s))
        or counts.sum() != len(s)
        or np.any(counts <= 0)
        or not np.isfinite(reduced).all()
        or not np.isfinite(s).all()
        or not np.isfinite(temperature)
        or temperature <= 0
    ):
        raise ValueError("Finite samples, positive counts and temperature required")
    if (
        edges.ndim != 1
        or len(edges) < 3
        or not np.isfinite(edges).all()
        or np.any(np.diff(edges) <= 0)
        or s.min() < edges[0]
        or s.max() >= edges[-1]
    ):
        raise ValueError("Increasing finite edges must strictly cover sample maxima")
    centers = (edges[:-1] + edges[1:]) / 2
    bins = np.searchsorted(edges, s, side="right") - 1
    nbin = len(centers)
    population = np.bincount(bins, minlength=nbin)
    occupied = population > 0
    reference = np.flatnonzero(np.isclose(centers, reference_s, rtol=0, atol=1e-12))
    if len(reference) != 1 or not occupied[reference[0]]:
        raise ValueError("Reference must identify one occupied histogram bin")
    fes = FES(reduced, counts, mbar_options={"solver_protocol": "robust"})
    fes.generate_fes(np.zeros(len(s)), s, histogram_parameters={"bin_edges": [edges]})
    estimate = fes.get_fes(
        centers[occupied],
        reference_point="from-specified",
        fes_reference=float(reference_s),
        uncertainty_method="analytical",
    )
    scale = KB * temperature
    values, errors = np.full(nbin, np.nan), np.full(nbin, np.nan)
    values[occupied] = estimate["f_i"] * scale
    # FES histogram estimates probability mass; this PMF reports density.
    values[occupied] += scale * np.log(
        np.diff(edges)[occupied] / np.diff(edges)[reference[0]]
    )
    errors[occupied] = estimate["df_i"] * scale
    if (
        not np.isfinite(values[occupied]).all()
        or not np.isfinite(errors[occupied]).all()
    ):
        raise ValueError("PyMBAR returned nonfinite estimates or uncertainties")
    mbar = fes.get_mbar()
    logweights = -logsumexp(
        np.log(counts)[:, None] + mbar.f_k[:, None] - reduced, axis=0
    )
    local_ess = np.zeros(nbin)
    for i in np.flatnonzero(occupied):
        local = logweights[bins == i]
        normalized = np.exp(local - logsumexp(local))
        local_ess[i] = 1 / np.square(normalized).sum()
    owners = np.repeat(np.arange(len(counts)), counts)
    contributors = np.zeros(nbin, int)
    contributors[np.unique(bins)] = [
        len(np.unique(owners[bins == i])) for i in np.unique(bins)
    ]

    def nullable(array: np.ndarray) -> list[float | None]:
        return [float(value) if np.isfinite(value) else None for value in array]

    return dict(
        centers=centers.tolist(),
        edges=edges.tolist(),
        reference_s=reference_s,
        pmf_kcal_mol=nullable(values),
        standard_error_kcal_mol=nullable(errors),
        lower_95_kcal_mol=nullable(values - 1.959963984540054 * errors),
        upper_95_kcal_mol=nullable(values + 1.959963984540054 * errors),
        selected_frames_per_bin=population.tolist(),
        weight_ess_per_bin=local_ess.tolist(),
        contributing_windows_per_bin=contributors.tolist(),
        confidence_method="PyMBAR histogram asymptotic covariance; decorrelated samples",
        interval_scope="Pointwise approximate normal 95%; conditional statistical error only",
        sampling_converged=False,
        global_barrier_accepted=False,
        mbar_f_k=mbar.f_k.tolist(),
        overlap=mbar.compute_overlap()["matrix"].tolist(),
        unbiased_ESS=float(
            np.exp(2 * logsumexp(logweights) - logsumexp(2 * logweights))
        ),
    )


def report_uncertainty(
    config: dict[str, Any],
    blocks: list[dict[str, Any]],
    report: dict[str, Any],
    *,
    bootstrap_replicates: int = 0,
    bootstrap_workers: int = 1,
) -> dict[str, Any]:
    """Add errors without changing sampling, bins, reference, or source estimator.

    Args:
        config: Frozen source analysis request.
        blocks: Hash-validated blocks accepted by the source analysis.
        report: Frozen source analysis result, which remains unchanged.
        bootstrap_replicates: Explicit completed-block refits; zero disables them.
        bootstrap_workers: Bounded CPU processes for optional bootstrap refits.

    Returns:
        Separate uncertainty/support report after estimator identity checks.
    """
    if config["temperature_K"] != report["temperature_K"] or report["pmf"] is None:
        raise ValueError("A matching physical temperature and source PMF are required")
    samples, counts = selected_samples(config["windows"], blocks, report)
    s = project_path(samples[:, :3], np.asarray(config["path_3cv_A"]))["s"]
    if not np.allclose(s, report["projection_diagnostics"]["s"], rtol=0, atol=1e-12):
        raise ValueError("Frozen source path projections changed")
    biases = [RelayBiasWindow(**w["bias"]) for w in config["windows"]]
    raw = dict(zip(("q_acyl", "q_act", "q_del"), samples[:, :3].T, strict=True))
    reduced = cross_bias_energies(raw, biases) / (KB * config["temperature_K"])
    reference_s = report["checkpoint"]["plot_common_es_reference_s"]
    result = histogram_uncertainty(
        reduced,
        counts,
        s,
        np.array(config["s_bin_edges"]),
        config["temperature_K"],
        reference_s,
    )
    original = np.array(report["pmf"]["pmf_kcal_mol"], float)
    centers = np.array(result["centers"])
    original -= original[np.flatnonzero(np.isclose(centers, reference_s))[0]]
    if not np.allclose(
        original, np.array(result["pmf_kcal_mol"], float), equal_nan=True, atol=1e-6
    ):
        raise ValueError("Uncertainty estimator failed frozen-PMF identity check")
    result["source_estimator_identity_verified"] = True
    result["temperature_K"] = config["temperature_K"]
    result["N_k"] = counts.tolist()
    if bootstrap_replicates:
        from .relay_pmf_block_bootstrap import block_bootstrap

        lengths = [
            [
                summary["selected_frames"]
                for summary in report["blocks"]
                if summary["window"] == window["label"]
            ]
            for window in config["windows"]
        ]
        result["block_bootstrap"] = block_bootstrap(
            reduced,
            counts,
            s,
            np.array(config["s_bin_edges"]),
            config["temperature_K"],
            reference_s,
            lengths,
            replicates=bootstrap_replicates,
            seed=20261126,
            workers=bootstrap_workers,
        )
    return result
