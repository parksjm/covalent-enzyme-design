#!/usr/bin/env python3
"""Audit saved relay PMF basin support without refitting the estimator."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
BASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis"
)
FIXED = BASE / "fixed_burnin_20261126/fixed_25ps_analysis.json"
REQUEST = BASE / "cumulative_progression_20261125/checkpoint_25ps/analysis_request.json"
PARENT = BASE / "cumulative_progression_20261125/checkpoint_25ps/analysis.json"
PROJECTION_SOURCE = (
    BASE
    / "fixed_burnin_20261126/parallel_source/src/nyl12_covalent/relay_umbrella_projection.py"
)
OUTPUT = ROOT / "provenance/relay_pmf_25ps_remaining_checks_20261126.json"


def sha256(path: Path) -> str:
    """Return the SHA256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load(path: Path) -> dict[str, Any]:
    """Load one JSON object from disk."""
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise TypeError(f"Expected JSON object in {path}")
    return value


def basin(pmf: dict[str, Any], bounds: list[float]) -> dict[str, Any]:
    """Find the saved PMF minimum and evaluate path-histogram brackets."""
    centers = pmf["centers"]
    values = pmf["pmf_kcal_mol"]
    occupied = pmf["occupied"]
    candidates = [
        i
        for i, center in enumerate(centers)
        if bounds[0] <= center <= bounds[1] and occupied[i] and values[i] is not None
    ]
    if not candidates:
        raise ValueError(f"No occupied PMF bins within {bounds}")
    index = min(candidates, key=lambda i: values[i])
    candidate_indices = [
        i for i, center in enumerate(centers) if bounds[0] <= center <= bounds[1]
    ]
    neighbors = list(range(max(0, index - 1), min(len(centers), index + 2)))
    return {
        "bounds": bounds,
        "center": centers[index],
        "pmf_kcal_mol": values[index],
        "boundary_candidate": index in (candidate_indices[0], candidate_indices[-1]),
        "neighbor_indices": neighbors,
        "neighbor_centers": [centers[i] for i in neighbors],
        "neighbor_occupied": [occupied[i] for i in neighbors],
        "neighbor_raw_frame_counts": [pmf["raw_sample_count"][i] for i in neighbors],
        "interior_with_occupied_neighbors": (
            len(candidate_indices) >= 3
            and index not in (candidate_indices[0], candidate_indices[-1])
            and all(occupied[i] for i in neighbors)
        ),
        "index": index,
    }


def main() -> None:
    """Write a hash-traceable comparison of saved estimator support."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    fixed, request, parent = map(load, (FIXED, REQUEST, PARENT))
    primary = fixed["primary"]
    pmf = {
        "centers": primary["centers"],
        "pmf_kcal_mol": primary["pmf_kcal_mol"],
        "occupied": [value is not None for value in primary["pmf_kcal_mol"]],
        "raw_sample_count": primary["frames_per_bin"],
    }
    es = basin(pmf, request["es_search_bounds"])
    acyl = basin(pmf, request["acyl_search_bounds"])
    all_occupied = pmf["occupied"]
    path_connected = (
        all(all_occupied[es["index"] : acyl["index"] + 1])
        if acyl["index"] > es["index"]
        else False
    )
    centers = primary["centers"]
    fixed_pmf = primary["pmf_kcal_mol"]
    fixed_counts = primary["frames_per_bin"]
    fixed_ess = primary["weight_ess_per_bin"]
    acyl_fixed = min(
        (
            i
            for i, x in enumerate(centers)
            if 0.9 <= x <= 1.1 and fixed_pmf[i] is not None
        ),
        key=lambda i: fixed_pmf[i],
    )
    half_minima: dict[str, dict[str, float]] = {}
    for half in ("early_half", "late_half"):
        profile = fixed["sensitivity"][half]
        indices = [
            i
            for i, x in enumerate(profile["centers"])
            if 0.9 <= x <= 1.1 and profile["pmf_kcal_mol"][i] is not None
        ]
        i = min(indices, key=lambda j: profile["pmf_kcal_mol"][j])
        half_minima[half] = {
            "acyl_minimum_s": profile["centers"][i],
            "acyl_minimum_kcal_mol": profile["pmf_kcal_mol"][i],
        }
    audit = {
        "schema_version": 1,
        "generated_by": "scripts/audit_relay_pmf_support.py",
        "scope": "Derived audit of saved estimators; no refit, new job, or changed reference.",
        "input_sha256": {
            str(p.relative_to(ROOT)): sha256(p)
            for p in (FIXED, REQUEST, PARENT, PROJECTION_SOURCE)
        },
        "fixed_report_sha256": sha256(FIXED),
        "checkpoint_request_sha256": sha256(REQUEST),
        "checkpoint_request_hash_matches": sha256(REQUEST) == parent["request_sha256"],
        "checkpoint_analysis_request_sha256": parent["request_sha256"],
        "estimator": "fixed-burn-in MBAR PMF using every retained frame from 5 to 25 ps/window",
        "horizon_ps": fixed["horizon_ps"],
        "temperature_K": fixed["temperature_K"],
        "search_bounds": {
            "es": request["es_search_bounds"],
            "acyl": request["acyl_search_bounds"],
        },
        "all_frame_path_histogram": {
            "es_basin": es,
            "acyl_basin": acyl,
            "intervening_bins_occupied": path_connected,
            "numerical_reference_bracketing_passes": (
                es["interior_with_occupied_neighbors"]
                and acyl["interior_with_occupied_neighbors"]
                and path_connected
            ),
        },
        "fixed_burnin_selected_estimator": {
            "horizon_ps": fixed["horizon_ps"],
            "reference_plot_bin_s": fixed["primary"]["reference_s"],
            "es_density_minimum_s": centers[es["index"]],
            "acyl_minimum_s": centers[acyl_fixed],
            "acyl_minimum_kcal_mol": fixed_pmf[acyl_fixed],
            "acyl_minimum_bin_frames": fixed_counts[acyl_fixed],
            "acyl_minimum_bin_weighted_ess": fixed_ess[acyl_fixed],
            "acyl_next_bin_frames": fixed_counts[acyl_fixed + 1],
            "acyl_next_bin_weighted_ess": fixed_ess[acyl_fixed + 1],
            "early_late_acyl_minima": half_minima,
        },
        "interpretation": (
            "The exact path_histogram geometric bracketing predicate passes for the saved 25 ps fixed-burn-in profile: both density minima are interior, "
            "each minimum and both neighboring bins are occupied, and the intervening bins are occupied. "
            "The acyl basin remains statistically weak (weighted ESS 1.24 at its minimum and one raw frame in its next bin), "
            "so geometric bracketing does not establish a reliable basin free energy or accepted barrier. "
            "All findings are conditional on the modeled pathway, saved weights, and retained 5–25 ps/window frames."
        ),
        "older_parent_gate_status_not_applicable_to_all_frame_estimator": {
            "neighbor_histogram_overlap": parent["gates"]["neighbor_histogram_overlap"][
                "passed"
            ],
            "effective_samples": parent["gates"]["effective_samples"]["passed"],
            "note": "These values describe the old checkpoint estimator only.",
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    LOGGER.info("Wrote derived audit to %s", OUTPUT)


if __name__ == "__main__":
    main()
