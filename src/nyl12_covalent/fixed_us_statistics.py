"""Shared fixed-window selection, native-bias MBAR, and block uncertainty."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

from .neb_umbrella_analysis import KB
from .relay_fixed_burnin import fixed_window_data
from .relay_fixed_pmf import contiguous_bootstrap, fit_samples
from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from .relay_umbrella_projection import project_path


def _features(config: dict[str, Any], values: np.ndarray) -> dict[str, np.ndarray]:
    """Name production columns and derive legacy q_h when explicitly biased."""
    names = tuple(
        config.get(
            "sample_feature_names",
            ("q_acyl", "q_act", "q_del", *config.get("diagnostic_columns", ())),
        )
    )
    if values.ndim != 2 or values.shape[1] != len(names):
        raise ValueError("Production columns do not match declared feature names")
    features = {name: values[:, i] for i, name in enumerate(names)}
    needed = {
        name
        for window in config["windows"]
        for name in window["bias"].get("feature_names", ("q_acyl", "q_h"))
    }
    if "q_h" in needed and "q_h" not in features:
        if not {"q_act", "q_del"} <= features.keys():
            raise ValueError("q_h bias requires q_h or both q_act and q_del")
        features["q_h"] = features["q_act"] + features["q_del"]
    return features


def fit_fixed_data(
    config: dict[str, Any],
    blocks: list[dict[str, Any]],
    burn_in_ps: float,
    end_ps: float,
) -> tuple[dict[str, Any], tuple, list[dict[str, Any]]]:
    """Select all frames after a fixed cut and fit the native-bias PMF.

    Args:
        config: Window biases, ordered route reference, bins and fixed reference.
        blocks: Accepted continuous production segments with timestamps.
        burn_in_ps: Per-window production burn-in boundary.
        end_ps: Inclusive per-window production horizon.

    Returns:
        Point-estimate report, immutable MBAR/bootstrap arguments, and selection
        diagnostics. A direct two-CV bias may use q_act + q_del as q_h while
        retaining the native feature names and centers.
    """
    values, counts, summaries = fixed_window_data(
        config["windows"], blocks, burn_in_ps=burn_in_ps, end_ps=end_ps
    )
    biases = [RelayBiasWindow(**window["bias"]) for window in config["windows"]]
    features = _features(config, values)
    reduced = cross_bias_energies(features, biases) / (KB * config["temperature_K"])
    coordinate = project_path(values[:, :3], np.asarray(config["path_3cv_A"]))["s"]
    arguments = (
        reduced,
        counts,
        coordinate,
        np.asarray(config["s_bin_edges"], dtype=float),
        config["temperature_K"],
        config["fixed_reference_s"],
    )
    report = fit_samples(*arguments)
    report.update(burn_in_ps=burn_in_ps, end_ps=end_ps)
    return report, arguments, summaries


def run_fixed_statistics(
    config: dict[str, Any],
    blocks: list[dict[str, Any]],
    horizon: float,
    replicates: int,
    workers: int,
) -> dict[str, Any]:
    """Run the shared primary fit, 5/2.5-ps bootstraps, and sensitivities.

    Args:
        config: Validated native biases, three-CV reference, bins and temperature.
        blocks: Accepted uninterrupted production segments for every window.
        horizon: Common production horizon in picoseconds, before the discard.
        replicates: Number of full MBAR refits for each bootstrap block length.
        workers: Explicit bootstrap worker count, between one and eight.

    Returns:
        Primary estimate, bootstrap intervals, sensitivities and frame selections.
    """
    logging.info("Fitting fixed-discard PMF through %g ps/window", horizon)
    primary, arguments, summaries = fit_fixed_data(config, blocks, 5.0, horizon)
    if not primary["connected"]:
        raise ValueError("Central MBAR network is disconnected; no global PMF")
    intervals = {round(row["frame_interval_ps"], 8) for row in summaries}
    if len(intervals) != 1:
        raise ValueError("Common sampling interval required for chunk resampling")
    interval = intervals.pop()
    bootstraps = {}
    for chunk in (5.0, 2.5):
        ratio = chunk / interval
        if not np.isclose(ratio, round(ratio), atol=1e-6):
            raise ValueError("Chunk size does not match sample spacing")
        logging.info("Starting %g-ps block bootstrap: %d refits", chunk, replicates)
        result = contiguous_bootstrap(
            *arguments,
            chunk_frames=round(ratio),
            replicates=replicates,
            workers=workers,
            seed=20261126,
        )
        bootstraps[f"{chunk:g}ps"] = result | {"chunk_ps": chunk}
        logging.info("Completed %g-ps block bootstrap", chunk)
    sensitivity = {}
    for cut in (0.0, 2.5, 10.0):
        fit, _, _ = fit_fixed_data(config, blocks, cut, horizon)
        sensitivity[f"burn_{cut:g}ps"] = fit
    midpoint = (5.0 + horizon) / 2
    for name, lower, upper in (
        ("early_half", 5.0, midpoint),
        ("late_half", midpoint, horizon),
    ):
        fit, _, _ = fit_fixed_data(config, blocks, lower, upper)
        sensitivity[name] = fit
    return dict(
        primary=primary,
        bootstrap=bootstraps,
        sensitivity=sensitivity,
        window_selection=summaries,
    )
