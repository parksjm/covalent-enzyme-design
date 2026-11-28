"""Diagnose SCC accuracy/step scaling, then conditionally complete validation."""

from __future__ import annotations

import argparse
import gzip
import logging
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

from nyl12_covalent.stationary_native import NativeState, write_json
from nyl12_covalent.stationary_precision import column_difference

LOGGER = logging.getLogger(__name__)


def main() -> None:
    """Measure derivative accuracy on one fixed chemical coordinate first."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    root = args.root.resolve()
    state = NativeState(root)
    column = 6562
    products = []
    errors = []
    records = []
    for step in (0.008, 0.004, 0.002, 0.001, 0.0005, 0.00025, 0.000125, 0.0000625):
        output = root / "accuracy_diagnostic" / f"step_{step}"
        product = column_difference(state, column, step, output)
        products.append(product)
        records.append(
            {"step_angstrom": step, "column_norm": float(np.linalg.norm(product))}
        )
        if len(products) > 1:
            error = float(np.linalg.norm(products[-1] - products[-2]))
            errors.append(error)
            records[-1].update(
                {
                    "absolute_difference_from_previous": error,
                    "relative_difference_from_previous": error
                    / np.linalg.norm(product),
                }
            )
        LOGGER.info("Completed chemical column at h=%s", step)
        write_json(
            root / "accuracy_diagnostic/progress.json", {"step_records": records}
        )
    ratio = errors[-1] / errors[-2] if errors[-2] else 0.0
    # A ratio near 1/4 supports central-difference truncation dominance.
    # This bounded screen cannot establish a global Hessian error bound.
    stable_indices = [
        index
        for index in range(2, len(records))
        if records[index]["relative_difference_from_previous"] < 1e-4
        and records[index - 1]["relative_difference_from_previous"] < 1e-4
    ]
    stable = bool(stable_indices)
    selected_index = (
        min(
            stable_indices,
            key=lambda index: records[index]["absolute_difference_from_previous"],
        )
        if stable
        else len(records) - 1
    )
    selected_step = records[max(1, selected_index - 1)]["step_angstrom"]

    electronic = []
    for path in sorted((root / "accuracy_diagnostic").glob("**/*.mdout.gz")):
        with gzip.open(path, "rt") as stream:
            text = stream.read()
        rows = re.findall(
            r"QMMM SCC-DFTB:\s+(\d+)\s+([-+\d.E]+)\s+([-+\d.E]+)\s+([-+\d.E]+)\s+\d+",
            text,
        )
        electronic.append(
            {
                "file": str(path.relative_to(root)),
                "last_verbose_scc_row": rows[-1] if rows else None,
                "failure_warning": "Convergence could not be achieved" in text,
            }
        )
    if any(
        row["last_verbose_scc_row"] is None or row["failure_warning"]
        for row in electronic
    ):
        stable = False
    result = {
        "precision_screen_passed": stable,
        "selected_hessian_step_angstrom": selected_step,
        "selection_basis": "Intermediate step in at least three-point stable derivative plateau; ratio alone is not required once numerical floor is reached.",
        "column": column,
        "step_records": records,
        "difference_ratio_finest_two": ratio,
        "scfconv": 1e-12,
        "tight_p_conv": 1,
        "density_conv": 1e-12,
        "original_density_conv": 5e-7,
        "electronic_records": electronic,
        "scope": "Same geometry,0K,DFTB3/MM,protonation,QM selection and4165movableatoms; electronic precision and derivative step only.",
        "interpretation": "Selected-column derivative scaling, not full-Hessian index or an operator-norm error bound.",
    }
    write_json(root / "accuracy_diagnostic/scientific_acceptance.json", result)
    if not stable:
        raise RuntimeError(
            "Measured derivative stability/SCC gate failed; full matrix not started"
        )
    entry = Path(__file__).with_name("run_stationary_validation.py")
    subprocess.run(
        [
            sys.executable,
            str(entry),
            "--root",
            str(root),
            "--workers",
            str(args.workers),
            "--step",
            str(selected_step),
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
