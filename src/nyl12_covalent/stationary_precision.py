"""Precision and serialization parity gates for native Hessian calculations."""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np

from nyl12_covalent.stationary_native import NativeState, write_json
from scripts.prepare_dimer_continuation import read_xyz_frames


def column_difference(
    state: NativeState, column: int, step: float, evidence: Path
) -> np.ndarray:
    """Return one mass-weighted central analytic-gradient difference."""
    _, plus = state.gradient(state.displaced(column, step), evidence / "plus")
    _, minus = state.gradient(state.displaced(column, -step), evidence / "minus")
    return ((plus - minus).reshape(-1) / (2 * step)) * (
        state.inv_sqrt * state.inv_sqrt[column]
    )


def precision_gate(root: Path, step: float) -> dict[str, object]:
    """Verify native old/new parity, force convergence and derivative stability.

    Args:
        root: Immutable validation campaign directory.
        step: Cartesian derivative step for the full Hessian.

    Returns:
        Actual numerical gate results. Failure blocks full Hessian calculation.
    """
    state = NativeState(root)
    out = root / "precision"
    out.mkdir()
    old_energy, old_gradient = state.gradient(
        state.xyz, out / "original", precision=False
    )
    new_energy, new_gradient = state.gradient(state.xyz, out / "precise")
    rounding_bound = 0.50001e-7 / 0.529177**2
    parity = float(np.max(np.abs(old_gradient - new_gradient)))
    mode = read_xyz_frames(root / "inputs/path_tsmode.xyz", len(state.indices))[-1]
    chemical_column = int(np.argmax(np.abs(mode.reshape(-1))))
    columns = sorted({chemical_column, 0, len(state.indices) * 3 - 1})
    differences = []
    relative_differences = []
    for column in columns:
        coarse = column_difference(state, column, step, out / f"column_{column}_coarse")
        fine = column_difference(state, column, step / 2, out / f"column_{column}_fine")
        differences.append(float(np.linalg.norm(coarse - fine)))
        relative_differences.append(
            differences[-1] / max(float(np.linalg.norm(fine)), 1e-14)
        )
    # A column probe is a precision diagnostic, not an operator-norm error bound.
    force = float(np.max(np.abs(new_gradient)) * 0.529177)
    reference = (root / "inputs/native_candidate.mdout").read_text()
    converged = bool(re.search(r"Final converged energy:", reference))
    results = {
        "passed": bool(
            parity <= rounding_bound
            and abs(old_energy - new_energy) <= 1e-9
            and force <= 1e-4
            and max(relative_differences) <= 1e-4
            and converged
        ),
        "active_atoms": int(state.active.sum()),
        "dimension": int(state.active.sum() * 3),
        "electronic_temperature_kelvin": 0,
        "qm_charge": -2,
        "energy_difference_hartree": new_energy - old_energy,
        "serialization_parity_max_gradient_difference": parity,
        "old_rounding_bound_hartree_per_angstrom": rounding_bound,
        "max_actual_gradient_hartree_per_bohr": force,
        "native_midpoint_converged": converged,
        "columns_checked": columns,
        "step_angstrom": step,
        "half_step_column_difference_norms": differences,
        "relative_half_step_column_differences": relative_differences,
        "engineering_relative_derivative_stability_gate": 1e-4,
        "gate_meaning": "Dimensionless derivative stability screen only; not a universal curvature tolerance or full-Hessian sign guarantee.",
        "limitation": "Three-column step check does not bound full Hessian truncation error.",
    }
    write_json(root / "precision_acceptance.json", results)
    return results
