"""Run a gated full-shell Nyl12 stationary-point validation loop."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from nyl12_covalent.stationary_endpoints import endpoints
from nyl12_covalent.stationary_hessian import full_hessian
from nyl12_covalent.stationary_native import write_json
from nyl12_covalent.stationary_precision import precision_gate


def main() -> None:
    """Execute parity, Hessian inertia, downhill and endpoint-minimum gates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--step", type=float, default=0.002)
    parser.add_argument("--precision-only", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not 1 <= args.workers <= 8 or not 0 < args.step <= 0.005:
        raise ValueError("Invalid authorized resource or derivative settings")
    root = args.root.resolve()
    path = root / "precision_acceptance.json"
    if path.exists():
        precision = json.loads(path.read_text())
    else:
        precision = precision_gate(root, args.step)
    if not precision["passed"]:
        raise RuntimeError("Precision/parity gate failed; full Hessian not launched")
    if args.precision_only:
        return
    result = full_hessian(
        root,
        root / "inputs/optimized.ncrst",
        root / "ts_hessian",
        args.workers,
        args.step,
        1,
    )
    if not result["passed"]:
        raise RuntimeError("Complete Hessian has not established numerical index one")
    receipt = json.loads((root / "precision/precise/receipt.json").read_text())
    restarts = endpoints(root, float(receipt["energy_hartree"]))
    minima = []
    for restart in restarts:
        minima.append(
            full_hessian(
                root,
                restart,
                restart.parent / "full_hessian",
                args.workers,
                args.step,
                0,
            )
        )
    write_json(
        root / "scientific_acceptance.json",
        {
            "passed": bool(all(item["passed"] for item in minima)),
            "saddle_index_one": result["passed"],
            "both_downhill_endpoints_connected": True,
            "endpoint_minimum_inertia": [item["passed"] for item in minima],
            "interpretation": "Same 0 K DFTB3/MM potential and full original 4165-atom movable region. Numerical Hessian inertia remains subject to recorded finite-difference diagnostics; no PMF or experimental mechanism claim.",
        },
    )


if __name__ == "__main__":
    main()
