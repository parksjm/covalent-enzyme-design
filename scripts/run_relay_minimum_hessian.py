"""Continue an accepted fine-step diagnostic into full-shell minimum inertia."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.relay_minimum import journal, verify_frozen
from nyl12_covalent.relay_minimum_run import CAVEAT
from nyl12_covalent.stationary_hessian import full_hessian
from nyl12_covalent.stationary_native import digest, write_json
from scripts.watch_relay_minimum import resource_gate

LOGGER = logging.getLogger(__name__)
STEP = 0.000125


def qualify(report: dict) -> None:
    """Reject a missing plateau, changed criterion or absent sampled step.

    Args:
        report: Retained fine-step diagnostic, not a minimum classification.
    """
    if (
        report["dimension"] != 12495
        or report["engineering_relative_column_gate"] != 1e-4
        or not report["two_adjacent_fine_pair_passes"]
        or not report["HB_yz_reciprocity_improved_at_least_tenfold"]
    ):
        raise ValueError("Fine-step precision qualification failed")
    pairs = report["adjacent_sampled_step_comparisons"]
    selected = next(item for item in pairs if item["coarse_step_A"] == STEP)
    if selected["fine_step_A"] != STEP / 2 or not selected["all_nine_pass"]:
        raise ValueError("Selected full-Hessian step is not on the tested plateau")


def main() -> None:
    """Validate immutable evidence, then compute the original full-shell Hessian."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not 1 <= args.workers <= 8:
        raise ValueError("Require one to eight single-CPU workers")
    root = args.root.resolve(strict=True)
    source = root / "05_fine_step_diagnostic/diagnostic_acceptance.json"
    report = json.loads(source.read_text())
    qualify(report)
    for relative, expected in report["source_column_sha256"].items():
        if digest(root / relative) != expected:
            raise ValueError(f"Fine-step evidence changed: {relative}")
    frozen = root / "01_frozen"
    verify_frozen(frozen)
    if not resource_gate(root, args.workers)["passed"]:
        raise RuntimeError("Full-Hessian resource headroom unavailable")
    LOGGER.info("Fine-step evidence accepted: %s; h=%g A", digest(source), STEP)
    if args.preflight_only:
        return
    summary = {
        "status": "full_hessian_running",
        "numerically_verified": False,
        "selected_step_A": STEP,
        "fine_step_report_sha256": digest(source),
        "selection": "Largest sampled step on the accepted fine-step plateau",
        "caveat": CAVEAT,
    }
    target = root / "verification_summary.json"
    write_json(target, summary)
    journal(
        "intent", "Run full12495-DOF minimum Hessian h=.000125A expected_negative=0"
    )
    try:
        result = full_hessian(
            frozen,
            frozen / "inputs/optimized.ncrst",
            root / "06_full_hessian",
            args.workers,
            STEP,
            expected_negative=0,
        )
    except Exception as error:
        summary.update(status="incomplete", error=str(error))
        write_json(target, summary)
        journal("outcome", f"Full-shell minimum verification incomplete: {error}")
        raise
    summary.update(
        status="numerical_local_minimum" if result["passed"] else "not_verified",
        numerically_verified=bool(result["passed"]),
        literal_inertia=result,
    )
    write_json(target, summary)
    journal(
        "outcome",
        f"Full-shell numerical minimum acceptance={result['passed']}",
        (target,),
    )


if __name__ == "__main__":
    main()
