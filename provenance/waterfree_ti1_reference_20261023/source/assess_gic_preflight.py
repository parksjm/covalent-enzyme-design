"""Check real-engine preservation of frozen Cartesian and torsion coordinates."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import GetDihedralDeg

from nyl12_covalent.gaussian_checks import final_orientation


def main() -> None:
    """Require successful GIC toy optimization without treating it as TI1 data."""
    root = Path(__file__).resolve().parent / "gic_preflight_01"
    text = (root / "optimization.log").read_text()
    elements, xyz = final_orientation(text)
    start = np.asarray([[0, 0, 0], [1.4, 0, 0], [-0.4, 0.9, 0], [1.8, 0.5, 0.8]])
    shifts = np.linalg.norm(xyz - start, axis=1)

    def chi(points: np.ndarray) -> float:
        conformer = Chem.Conformer(4)
        for i, point in enumerate(points):
            conformer.SetAtomPosition(i, tuple(float(v) for v in point))
        return GetDihedralDeg(conformer, 2, 0, 1, 3)

    error = abs((chi(xyz) - chi(start) + 180) % 360 - 180)
    checks = {
        "normal_termination": "Normal termination of Gaussian 16" in text,
        "optimization_completed": "Optimization completed." in text,
        "atom_order": elements == [8, 8, 1, 1],
        "frozen_Cartesians": bool(max(shifts[:3]) <= 1e-6),
        "frozen_torsion": error <= 0.001,
        "mobile_atom_relaxed": bool(shifts[3] > 1e-4),
    }
    report = {
        "accepted": all(checks.values()),
        "checks": checks,
        "atom_shifts_A": shifts.tolist(),
        "torsion_error_deg": error,
        "scope": "HF/STO-3G installation/optimizer test only; not TI1 reference data",
    }
    (root / "scientific_acceptance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    logging.info("GIC frozen XYZ/torsion test accepted: %s", report["accepted"])
    if not report["accepted"]:
        raise SystemExit(3)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
