"""Explicit native optimizer displacement conversion and atom-step cap."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .native_sander_cli import ANGSTROM_PER_BOHR


def capped_step(step_bohr: np.ndarray, cap_A: float) -> np.ndarray:
    """Convert a finite native displacement to Å and cap each atom collectively.

    Args:
        step_bohr: Flat or (N,3) displacement supplied by native DL-FIND.
        cap_A: Positive maximum atom displacement in angstrom.

    Returns:
        Flat displacement in angstrom scaled uniformly to preserve direction.

    Raises:
        ValueError: If native displacement or cap is invalid.
    """
    step = np.asarray(step_bohr, dtype=float).reshape(-1, 3) * ANGSTROM_PER_BOHR
    if not np.isfinite(step).all() or not np.isfinite(cap_A) or cap_A <= 0:
        raise ValueError("Native step and positive atom cap must be finite")
    largest = float(np.linalg.norm(step, axis=1).max(initial=0))
    return (step * min(1.0, cap_A / largest) if largest else step).ravel()


def runtime_environment(amber: Path) -> dict[str, str]:
    """Return the explicit single-threaded environment for an Amber prefix."""
    return {
        "AMBERHOME": str(amber),
        "LD_LIBRARY_PATH": str(amber / "lib"),
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
