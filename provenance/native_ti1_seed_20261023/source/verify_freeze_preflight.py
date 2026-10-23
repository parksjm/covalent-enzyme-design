"""Verify the Gaussian partial-optimization syntax against a real water job."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

from nyl12_covalent.gaussian_checks import final_orientation


def main() -> None:
    """Require two frozen atoms to stay fixed while the third atom relaxes."""
    root = Path(__file__).resolve().parent / "freeze_flags_preflight_01"
    text = (root / "optimization.log").read_text()
    numbers, xyz = final_orientation(text)
    initial = np.array([[1, 2, 3], [1, 2.757, 3.586], [1, 1.243, 3.586]])
    shift = np.linalg.norm(xyz - initial, axis=1)
    checks = {
        "atom_order": numbers == [8, 1, 1],
        "normal_termination": "Normal termination of Gaussian 16" in text,
        "optimized": "Optimization completed." in text,
        "frozen_atoms": bool(max(shift[:2]) <= 1e-6),
        "mobile_atom_changed": bool(shift[2] > 1e-4),
    }
    report = {
        "accepted": all(checks.values()),
        "checks": checks,
        "atom_displacements_A": shift.tolist(),
        "scope": "Installation/freeze syntax only; HF/STO-3G water is not TI1 parameter data",
    }
    (root / "scientific_acceptance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    logging.info("Freeze syntax accepted: %s", report["accepted"])
    if not report["accepted"]:
        raise SystemExit(3)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
