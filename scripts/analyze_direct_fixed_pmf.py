#!/usr/bin/env python3
"""Reanalyze the NEB-seeded direct-US campaign with shared PMF statistics."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

import numpy as np

from nyl12_covalent.direct_us_fixed_io import DIAGNOSTICS, load_direct_dataset
from nyl12_covalent.fixed_us_statistics import run_fixed_statistics
from nyl12_covalent.relay_umbrella_analysis_io import digest


def main() -> None:
    """Validate native sources, refit MBAR and preserve route-specific diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("campaign", "augmented", "relay-request", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--horizon", type=float, required=True)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--continuation-plan", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve existing PMF evidence")
    config, blocks, sources = load_direct_dataset(
        args.campaign.resolve(),
        args.augmented.resolve(),
        args.relay_request.resolve(),
        continuation_plan=(
            args.continuation_plan.resolve() if args.continuation_plan else None
        ),
        horizon_ps=args.horizon,
    )
    statistics = run_fixed_statistics(
        config, blocks, args.horizon, replicates=args.replicates, workers=args.workers
    )
    chemistry = []
    for window, summary in zip(
        config["windows"], statistics["window_selection"], strict=True
    ):
        parts = sorted(
            [block for block in blocks if block["window"] == window["label"]],
            key=lambda block: float(block["original_rows"][0]["time_ps"]),
        )
        rows = [row for part in parts for row in part["original_rows"]]
        selected = [rows[index] for index in summary["selected_indices"]]
        owners = [int(row["hg1_owner"]) for row in selected]
        chemistry.append(
            dict(
                window=window["label"],
                frames=len(owners),
                hg1_nearest_nalpha_frames=owners.count(13905),
                hg1_nearest_og1_frames=owners.count(13916),
                hg1_nearest_leaving_n_frames=owners.count(19929),
                frames_censored=0,
            )
        )
    result = dict(
        **statistics,
        protocol="fixed5ps_once_per_continuous_production_window_all_frames",
        horizon_ps=args.horizon,
        windows=len(config["windows"]),
        temperature_K=config["temperature_K"],
        native_heating_ps=1,
        native_equilibration_ps=5,
        additional_production_discard_ps=5,
        native_scc_policy={
            "legacy_scfconv": 1e-10,
            "full_step_scc_success_verified": False,
            "note": "Acceptance records reject reported SCC failures; the original stage protocol does not record full-step SCC verification.",
        },
        source_sha256=sources,
        config=config,
        campaign="Mechanism A: direct-NEB-seeded sampling; actual proton channels uncensored",
        projection=dict(
            features=["q_acyl", "q_act", "q_del"],
            metric="Equal-Angstrom Euclidean",
            path="Original24NEB images; continuous arclength; extended endpoint tangents",
            same_projection_algorithm_as_relay=True,
            route_specific_reference_path=True,
            three_cv_path_sha256=hashlib.sha256(
                np.asarray(config["path_3cv_A"], dtype="<f8").tobytes()
            ).hexdigest(),
        ),
        diagnostic_columns=["q_acyl", "q_act", "q_del", *DIAGNOSTICS],
        chemical_sampling=chemistry,
        source_code_sha256={
            str(path): digest(path)
            for path in (
                Path(__file__).resolve(),
                Path(__file__).resolve().parents[1]
                / "src/nyl12_covalent/direct_us_fixed_io.py",
                Path(__file__).resolve().parents[1]
                / "src/nyl12_covalent/direct_us_continuation_io.py",
                Path(__file__).resolve().parents[1]
                / "src/nyl12_covalent/fixed_us_statistics.py",
            )
        },
        sampling_converged=False,
        global_barrier_accepted=False,
        caveat="Matched statistics do not establish matched conformational ensembles or pure direct chemistry; both profiles remain conditional on sampled states.",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    logging.info("Wrote direct-US fixed-protocol PMF: %s", args.output)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
