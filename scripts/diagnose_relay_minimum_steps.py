"""Probe four finer derivative steps at the unchanged frozen relay structure."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.relay_minimum import journal, verify_frozen
from nyl12_covalent.relay_minimum_run import column_metrics, verify_prerequisites
from nyl12_covalent.stationary_hessian import worker_column
from nyl12_covalent.stationary_native import digest, write_json
from scripts.watch_relay_minimum import resource_gate

LOGGER = logging.getLogger(__name__)
COLUMNS = (0, 1, 2, 6585, 6586, 6587, 10449, 10450, 10451)
FINE_STEPS = (0.000125, 0.0000625, 0.00003125, 0.000015625)
LABELS = ("h_0p000125A", "h_0p0000625A", "h_0p00003125A", "h_0p000015625A")


def analyze(arrays: list[np.ndarray], steps: list[float]) -> dict[str, Any]:
    """Compare full column stability and sampled mixed-derivative reciprocity.

    Args:
        arrays: One 12495-by-9 mass-weighted column block per sampled step.
        steps: Corresponding decreasing displacements in Angstrom.

    Returns:
        Literal comparisons, unchanged stability gate and bounded interpretation.
    """
    if len(arrays) != len(steps) or len(arrays) < 3:
        raise ValueError("Require matched arrays and at least three steps")
    if any(
        a.shape != (12495, len(COLUMNS)) or not np.isfinite(a).all() for a in arrays
    ):
        raise ValueError("Require all nine finite complete active-space columns")
    if not all(np.isfinite(h) and h > 0 for h in steps) or any(
        a <= b for a, b in zip(steps, steps[1:])
    ):
        raise ValueError("Require strictly decreasing positive finite steps")
    comparisons = []
    reciprocity = []
    y, z = COLUMNS.index(6586), COLUMNS.index(6587)
    for coarse, fine, h, smaller in zip(arrays, arrays[1:], steps, steps[1:]):
        metrics = {
            str(c): column_metrics(coarse[:, i], fine[:, i])
            for i, c in enumerate(COLUMNS)
        }
        comparisons.append(
            {
                "coarse_step_A": h,
                "fine_step_A": smaller,
                "columns": metrics,
                "all_nine_pass": all(m["passed"] for m in metrics.values()),
            }
        )
    for block, step in zip(arrays, steps):
        yz, zy = float(block[6587, y]), float(block[6586, z])
        denominator = max(abs(yz), abs(zy))
        sampled = block[np.asarray(COLUMNS)]
        discrepancy = np.abs(sampled - sampled.T)
        reciprocity.append(
            {
                "step_A": step,
                "HB_Hzy": yz,
                "HB_Hyz": zy,
                "HB_yz_absolute_difference": abs(yz - zy),
                "HB_yz_relative_difference": abs(yz - zy) / denominator
                if denominator
                else None,
                "HB_yz_max_sampled_reciprocity_error": float(
                    discrepancy[[y, z], :].max()
                ),
            }
        )
    fine_pairs = [
        p["all_nine_pass"] for p in comparisons if p["coarse_step_A"] <= FINE_STEPS[0]
    ]
    plateau = any(a and b for a, b in zip(fine_pairs, fine_pairs[1:]))
    old, new = (
        reciprocity[0]["HB_yz_absolute_difference"],
        reciprocity[-1]["HB_yz_absolute_difference"],
    )
    improved = old > 0 and new <= old / 10
    return {
        "numerically_verified_minimum": False,
        "dimension": 12495,
        "column_units": "Hartree/(Angstrom^2 amu)",
        "engineering_relative_column_gate": 1e-4,
        "adjacent_sampled_step_comparisons": comparisons,
        "reciprocity_by_step": reciprocity,
        "two_adjacent_fine_pair_passes": plateau,
        "HB_yz_reciprocity_improved_at_least_tenfold": improved,
        "status": "ready_for_lead_review"
        if plateau and improved
        else "not_stable_or_reciprocity_unresolved",
        "limitation": "Selected columns only, not a full Hessian or minimum proof; no automatic matrix launch.",
    }


def run(root: Path, workers: int) -> dict[str, Any]:
    """Compute exactly36 new columns and analyze72 native gradient evaluations.

    Args:
        root: Existing w033 minimum_verification root with immutable01_frozen.
        workers: Positive maximum eight nice single-CPU native workers.

    Returns:
        Retained finite-step engineering diagnostics, never a minimum claim.
    """
    if not 1 <= workers <= 8:
        raise ValueError("Require1..8 single-CPU workers")
    stage = root / "05_fine_step_diagnostic"
    frozen = root / "01_frozen"
    journal(
        "intent",
        "Verify unchanged frozen/40prerequisite hashes and resources for72-gradient fine-step diagnostic",
    )
    verify_frozen(frozen)
    verify_prerequisites(root)
    resources = resource_gate(stage, workers)
    write_json(stage / "resource_preflight.json", resources)
    if not resources["passed"]:
        raise RuntimeError("Insufficient resources for bounded diagnostic")
    coarse = json.loads(
        (root / "03_step_probes/scientific_acceptance.json").read_text()
    )
    if set(coarse["columns"]) != {str(c) for c in COLUMNS}:
        raise ValueError("Original selected column identities changed")
    restart = frozen / "inputs/optimized.ncrst"
    write_json(
        stage / "specification.json",
        {
            "steps_A": list(FINE_STEPS),
            "columns_zero_based": list(COLUMNS),
            "native_gradient_calls": 72,
            "workers": workers,
            "restart_sha256": digest(restart),
            "library_sha256": digest(frozen / "library_build/libdlfind.so"),
            "relative_column_gate": 1e-4,
            "matrix_units": "Hartree/(Angstrom^2 amu)",
        },
    )
    journal(
        "outcome",
        "Exact frozen inputs and40prerequisite hashes retained; resource gate passed",
    )
    for label in LABELS:
        (stage / label).mkdir()
    tasks = [
        (frozen, stage / label, restart, c, h)
        for label, h in zip(LABELS, FINE_STEPS)
        for c in COLUMNS
    ]
    journal(
        "intent",
        "Evaluate36 full-length central columns at four finer steps;72 native analytic gradients",
    )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(worker_column, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), 1):
            future.result()
            write_json(
                stage / "progress.json",
                {"completed_columns": completed, "total_columns": len(tasks)},
            )
            LOGGER.info("Completed fine-step column%s/%s", completed, len(tasks))
    locations = [
        root / "03_step_probes/h_0p001A",
        root / "03_step_probes/h_0p0005A",
    ] + [stage / label for label in LABELS]
    arrays = [
        np.column_stack(
            [
                np.load(location / "columns" / f"column_{c:05d}" / "column.npy")
                for c in COLUMNS
            ]
        )
        for location in locations
    ]
    result = analyze(arrays, [0.001, 0.0005, *FINE_STEPS])
    result["source_column_sha256"] = {
        str(
            location.relative_to(root) / "columns" / f"column_{c:05d}" / "column.npy"
        ): digest(location / "columns" / f"column_{c:05d}" / "column.npy")
        for location in locations
        for c in COLUMNS
    }
    write_json(stage / "diagnostic_acceptance.json", result)
    journal(
        "outcome",
        f"Bounded fine-step diagnostic complete: {result['status']}; no full Hessian launched",
        (stage / "diagnostic_acceptance.json",),
    )
    return result


def main() -> None:
    """Run the authorized bounded diagnostic and retain failures explicitly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    try:
        run(args.root, args.workers)
    except Exception as error:
        journal(
            "outcome", f"Fine-step diagnostic failed: {type(error).__name__}: {error}"
        )
        write_json(
            args.root / "05_fine_step_diagnostic/failure.json",
            {"numerically_verified_minimum": False, "error": str(error)},
        )
        raise


if __name__ == "__main__":
    main()
