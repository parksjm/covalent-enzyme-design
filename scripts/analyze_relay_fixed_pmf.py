#!/usr/bin/env python3
"""Reanalyze a frozen relay checkpoint with fixed burn-in and block uncertainty."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np

from nyl12_covalent.relay_fixed_burnin import fixed_window_data
from nyl12_covalent.relay_fixed_pmf import contiguous_bootstrap, fit_samples
from nyl12_covalent.relay_umbrella_analysis import KB
from nyl12_covalent.relay_umbrella_analysis_io import digest, load_block
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from nyl12_covalent.relay_umbrella_projection import project_path


def prepare_fit(
    config: dict, blocks: list[dict], burn_in_ps: float, end_ps: float
) -> tuple[dict, tuple, list[dict]]:
    """Assemble unchanged physical-bias reweighting inputs after a fixed cut.

    Args:
        config: Hash-validated source request.
        blocks: All accepted source segments, with original timestamps.
        burn_in_ps: Time discarded once at the start of production per window.
        end_ps: Common production horizon, measured from production start.

    Returns:
        Point estimate, arguments for full MBAR refits, and selection summaries.
    """
    values, counts, summaries = fixed_window_data(
        config["windows"], blocks, burn_in_ps=burn_in_ps, end_ps=end_ps
    )
    biases = [RelayBiasWindow(**window["bias"]) for window in config["windows"]]
    features = dict(zip(("q_acyl", "q_act", "q_del"), values[:, :3].T, strict=True))
    reduced = cross_bias_energies(features, biases) / (KB * config["temperature_K"])
    projection = project_path(values[:, :3], np.array(config["path_3cv_A"]))
    arguments = (
        reduced,
        counts,
        projection["s"],
        np.array(config["s_bin_edges"]),
        config["temperature_K"],
        config["fixed_reference_s"],
    )
    report = fit_samples(*arguments)
    report.update(burn_in_ps=burn_in_ps, end_ps=end_ps)
    return report, arguments, summaries


def main() -> None:
    """Verify immutable sources and retain every fit, diagnostic and refit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replicates", required=True, type=int)
    parser.add_argument("--workers", required=True, type=int)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Retain prior evidence; choose a new output file")
    request_path = args.checkpoint / "analysis_request.json"
    parent_path = args.checkpoint / "analysis.json"
    config, parent = [
        json.loads(path.read_text()) for path in (request_path, parent_path)
    ]
    if digest(request_path) != parent["request_sha256"]:
        raise ValueError("Parent and request hash mismatch")
    blocks = [
        load_block(
            (request_path.parent / filename).resolve(strict=True),
            temperature=config["temperature_K"],
            hamiltonian=config["physical_hamiltonian_sha256"],
            diagnostic_columns=config["diagnostic_columns"],
        )
        for filename in config["accepted_block_receipts"]
    ]
    actual = [
        dict(
            receipt_sha256=block["receipt_sha256"], source_sha256=block["source_sha256"]
        )
        for block in blocks
    ]
    if actual != parent["accepted_sources"]:
        raise ValueError("Hash-validated blocks differ from accepted parent data")
    horizon = float(parent["checkpoint"]["actual_horizon_ps"])
    if horizon < 20:
        raise ValueError("This prespecified analysis requires at least20ps/window")
    config["fixed_reference_s"] = parent["checkpoint"]["plot_common_es_reference_s"]
    report, arguments, summaries = prepare_fit(config, blocks, 5.0, horizon)
    if not report["connected"]:
        raise ValueError("Central MBAR network is disconnected; no global PMF")
    dt = {round(row["frame_interval_ps"], 8) for row in summaries}
    if len(dt) != 1:
        raise ValueError("Common sampling interval required for chunk resampling")
    interval = dt.pop()
    bootstraps = {}
    for chunk in (5.0, 2.5):
        logging.info("Starting %gps block bootstrap: %s refits", chunk, args.replicates)
        ratio = chunk / interval
        if not np.isclose(ratio, round(ratio), atol=1e-6):
            raise ValueError("Chunk size does not match sample spacing")
        result = contiguous_bootstrap(
            *arguments,
            chunk_frames=round(ratio),
            replicates=args.replicates,
            workers=args.workers,
            seed=20261126,
        )
        bootstraps[f"{chunk:g}ps"] = result | dict(chunk_ps=chunk)
    sensitivities = {}
    for cut in (0.0, 2.5, 10.0):
        fit, _, _ = prepare_fit(config, blocks, cut, horizon)
        sensitivities[f"burn_{cut:g}ps"] = fit
    midpoint = (5.0 + horizon) / 2
    for name, lower, upper in (
        ("early_half", 5.0, midpoint),
        ("late_half", midpoint, horizon),
    ):
        fit, _, _ = prepare_fit(config, blocks, lower, upper)
        sensitivities[name] = fit
    files = [request_path, parent_path, Path(__file__)] + [
        Path(__file__).resolve().parents[1] / "src/nyl12_covalent" / filename
        for filename in ("relay_fixed_burnin.py", "relay_fixed_pmf.py")
    ]
    result = dict(
        protocol="fixed5ps_once_per_continuous_production_window_all_frames",
        horizon_ps=horizon,
        windows=len(config["windows"]),
        temperature_K=config["temperature_K"],
        native_heating_ps=1,
        native_equilibration_ps=5,
        additional_production_discard_ps=5,
        physical_hamiltonian_sha256=config["physical_hamiltonian_sha256"],
        path_3cv_sha256=parent["path_3cv_sha256"],
        source_sha256={str(path.resolve()): digest(path) for path in files},
        checkpoint=str(args.checkpoint.resolve()),
        accepted_sources=actual,
        window_selection=summaries,
        diagnostic_columns=["q_acyl", "q_act", "q_del"] + config["diagnostic_columns"],
        primary=report,
        bootstrap=bootstraps,
        sensitivity=sensitivities,
        caveat=(
            "Burn-in is a transparent provisional analysis choice, not proven equilibration. "
            "All-frame MBAR point estimate assumes stationary ergodic windows. "
            "Bootstrap intervals condition on visited configurations and block length; "
            "few blocks cannot exclude slower modes or initialization bias. "
            "No chemical censoring or changes to simulations, biases, path, bins or reference."
        ),
        sampling_converged=False,
        global_barrier_accepted=False,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    logging.info("Wrote fixed-burn-in PMF and sensitivity report: %s", args.output)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
