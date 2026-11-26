#!/usr/bin/env python3
"""Reanalyze a frozen relay checkpoint with fixed burn-in and block uncertainty."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.fixed_us_statistics import fit_fixed_data, run_fixed_statistics
from nyl12_covalent.relay_umbrella_analysis_io import digest, load_block


def prepare_fit(
    config: dict, blocks: list[dict], burn_in_ps: float, end_ps: float
) -> tuple[dict, tuple, list[dict]]:
    """Compatibility wrapper for shared fixed-window statistical preparation.

    Args:
        config: Hash-validated source request.
        blocks: All accepted source segments, with original timestamps.
        burn_in_ps: Time discarded once at the start of production per window.
        end_ps: Common production horizon, measured from production start.

    Returns:
        Point estimate, arguments for full MBAR refits, and selection summaries.
    """
    return fit_fixed_data(config, blocks, burn_in_ps, end_ps)


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
    statistics = run_fixed_statistics(
        config, blocks, horizon, args.replicates, args.workers
    )
    logging.info("Completed fixed-window block bootstraps: %s refits", args.replicates)
    files = [request_path, parent_path, Path(__file__)] + [
        Path(__file__).resolve().parents[1] / "src/nyl12_covalent" / filename
        for filename in (
            "relay_fixed_burnin.py",
            "relay_fixed_pmf.py",
            "fixed_us_statistics.py",
        )
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
        window_selection=statistics["window_selection"],
        diagnostic_columns=["q_acyl", "q_act", "q_del"] + config["diagnostic_columns"],
        primary=statistics["primary"],
        bootstrap=statistics["bootstrap"],
        sensitivity=statistics["sensitivity"],
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
