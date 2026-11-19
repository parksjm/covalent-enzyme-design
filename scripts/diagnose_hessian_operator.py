"""Check native DFTB3/MM Hessian-vector symmetry, step size and repeatability."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import re

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from nyl12_covalent.native_hessian_probe import finite_difference_hvp
from nyl12_covalent.native_sander_cli import (
    HARTREE_TO_KCAL_PER_MOL,
    NativeSanderEvaluator,
    load_calibration,
    prepare_native_input,
)
from scripts.prepare_dimer_continuation import read_xyz_frames

LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Hash one immutable diagnostic input.

    Args:
        path: Existing input file.

    Returns:
        Hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    """Run five bounded HVPs and save scalar diagnostic metrics."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("topology", "restart", "mdin", "mode", "calibration", "sander"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.output.exists():
        raise FileExistsError(args.output)
    for name in ("topology", "restart", "mdin", "mode", "calibration", "sander"):
        if not getattr(args, name).is_file():
            raise FileNotFoundError(getattr(args, name))
    scale, calibration = load_calibration(args.calibration)
    topology = AmberParm(str(args.topology))
    restart = Rst7.open(str(args.restart))
    coordinates = np.asarray(restart.coordinates, dtype=float).reshape((-1, 3))
    if len(coordinates) != len(topology.atoms):
        raise ValueError("Topology/restart atom counts disagree")
    mdin_source = args.mdin.read_text()
    match = re.search(r"\bactive\s*=\s*'([^']+)'", mdin_source, re.IGNORECASE)
    if match is None:
        raise ValueError("Missing active atom mask")
    active = np.asarray(AmberMask(topology, match.group(1)).Selection(), dtype=bool)
    indices = np.flatnonzero(active)
    if len(indices) != 4165:
        raise ValueError(f"Unexpected active atom count {len(indices)}")
    masses = np.asarray([atom.mass for atom in topology.atoms], dtype=float)[active]
    if np.any(masses <= 0) or not np.isfinite(masses).all():
        raise ValueError("Invalid active-atom masses")
    inverse_sqrt = 1.0 / np.sqrt(masses)
    mode = read_xyz_frames(args.mode, len(indices))[-1]
    v = (np.sqrt(masses[:, None]) * mode).reshape(-1)
    v /= np.linalg.norm(v)
    rng = np.random.default_rng(20261116)
    w = rng.normal(size=len(v))
    w -= np.dot(w, v) * v
    w /= np.linalg.norm(w)
    args.output.mkdir(parents=True)
    evaluation_dir = args.output / "evaluations"
    evaluation_dir.mkdir()
    evaluator = NativeSanderEvaluator(
        args.topology.resolve(),
        restart,
        active,
        prepare_native_input(mdin_source),
        args.sander.resolve(),
        evaluation_dir,
        scale,
    )

    def gradient(xyz: np.ndarray) -> np.ndarray:
        """Return calibrated native energy gradient in Hartree/Å."""
        _, kcal_gradient = evaluator(xyz)
        return kcal_gradient / HARTREE_TO_KCAL_PER_MOL

    def apply(direction: np.ndarray, step: float) -> np.ndarray:
        """Evaluate one central finite-difference mass-weighted HVP."""
        return finite_difference_hvp(
            direction, coordinates, indices, inverse_sqrt, gradient, step
        )

    av_005 = apply(v, 0.05)
    LOGGER.info("Completed H(v) at 0.05 Å√amu")
    aw_005 = apply(w, 0.05)
    av_010 = apply(v, 0.10)
    aw_010 = apply(w, 0.10)
    av_repeat = apply(v, 0.05)
    left, right = float(np.dot(v, aw_005)), float(np.dot(w, av_005))
    results = {
        "state": "0 K DFTB3/MM direct-route dimer candidate",
        "active_atoms": len(indices),
        "native_gradient_evaluations": len(evaluator.evidence),
        "hvp_steps_angstrom_sqrt_amu": [0.05, 0.10],
        "bilinear_symmetry": {
            "v_dot_Hw": left,
            "w_dot_Hv": right,
            "absolute_difference": abs(left - right),
            "relative_difference": abs(left - right)
            / max(abs(left), abs(right), 1.0e-12),
        },
        "step_sensitivity": {
            "chemical_mode_relative_norm": float(
                np.linalg.norm(av_005 - av_010) / np.linalg.norm(av_005)
            ),
            "random_mode_relative_norm": float(
                np.linalg.norm(aw_005 - aw_010) / np.linalg.norm(aw_005)
            ),
        },
        "repeatability": {
            "chemical_mode_relative_norm": float(
                np.linalg.norm(av_005 - av_repeat) / np.linalg.norm(av_005)
            )
        },
        "rayleigh_quotient_hartree_per_angstrom2_amu": {
            "chemical_mode": float(np.dot(v, av_005)),
            "random_mode": float(np.dot(w, aw_005)),
        },
        "force_calibration": calibration,
        "inputs_sha256": {
            str(getattr(args, name).resolve()): sha256(getattr(args, name))
            for name in ("topology", "restart", "mdin", "mode", "calibration", "sander")
        },
        "interpretation": "Operator check only; does not verify Hessian index.",
    }
    (args.output / "diagnosis.json").write_text(json.dumps(results, indent=2) + "\n")
    LOGGER.info("Completed %d native gradient evaluations", len(evaluator.evidence))


if __name__ == "__main__":
    main()
