"""Conditional dimer TS candidate from an accepted N-terminal relay profile."""

from __future__ import annotations

import math
import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np

from scripts.prepare_dimer_continuation import write_restart
from .relay_minimum import journal
from .stationary_native import NativeState, write_json


def _peak(
    rows: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]] | None:
    if (
        len(rows) < 3
        or rows[0].get("accepted") is not True
        or rows[-1].get("accepted") is not True
    ):
        raise ValueError("Profile requires accepted baseline and released endpoint")
    accepted = [row for row in rows if row.get("accepted") is True]
    energies = [float(row["physical_kcal_mol"]) for row in accepted]
    if not all(math.isfinite(value) for value in energies):
        raise ValueError("Nonfinite accepted physical energy")
    if len(accepted) < 3:
        return None
    index = max(range(1, len(accepted) - 1), key=lambda i: energies[i])
    if not (
        energies[index] > energies[index - 1] and energies[index] > energies[index + 1]
    ):
        return None
    return accepted[index - 1], accepted[index], accepted[index + 1]


def _option(source: str, name: str) -> str:
    values = re.findall(rf"\b{name}\s*=\s*'([^']+)'", source, re.IGNORECASE)
    if len(values) != 1:
        raise ValueError(f"Require exactly one {name} mask")
    return values[0]


def _dimer_input(source: str, reference: Path, mode: Path) -> str:
    parts = re.split(r"&dlfind\b", source, flags=re.IGNORECASE)
    if len(parts) != 2:
        raise ValueError("Require exactly one DL-FIND block")
    prefix = parts[0]
    if re.search(r"&wt\b|\bDISANG\s*=", prefix, re.IGNORECASE):
        raise ValueError("Reaction-coordinate restraints are not permitted")
    for key in ("ntr", "nmropt"):
        values = re.findall(rf"\b{key}\s*=\s*(\d+)", prefix, re.IGNORECASE)
        if values != ["0"]:
            raise ValueError(f"{key} must explicitly be zero")
    active, weights = _option(parts[1], "active"), _option(parts[1], "wtmask")
    qmmask = _option(prefix, "qmmask")
    if weights != qmmask:
        raise ValueError("Dimer weights differ from original QM mask")
    if "'" in str(reference) or "'" in str(mode):
        raise ValueError("Paths cannot contain quotes")
    return prefix + (
        "&dlfind\n dimer=.true., dimer_mode=1, dimer_tolrot=2., dimer_maxrot=20,\n"
        " optalg='LBFGS', crdrep='CART',\n"
        f" refcrd='{reference}', coords2(1)='{mode}',\n"
        f" wtmask='{weights}', active='{active}',\n"
        " tol=1.d-4, tole=1.d-7, maxcycle=5000, maxstep=0.05, verbosity=2,\n/\n"
    )


def _mode(
    previous: np.ndarray,
    peak: np.ndarray,
    following: np.ndarray,
    active: np.ndarray,
    reference: np.ndarray,
) -> np.ndarray:
    if (
        previous.shape != peak.shape
        or following.shape != peak.shape
        or reference.shape != peak.shape
        or peak.shape != (len(active), 3)
    ):
        raise ValueError("Restart coordinate inventories differ")
    if not all(np.isfinite(x).all() for x in (previous, peak, following, reference)):
        raise ValueError("Nonfinite restart coordinates")
    if any(
        not np.array_equal(x[~active], reference[~active])
        for x in (previous, peak, following)
    ):
        raise ValueError("A path restart changed the frozen exterior")
    tangent = following[active] - previous[active]
    length = float(np.linalg.norm(tangent))
    if not math.isfinite(length) or length <= 1e-12:
        raise ValueError("Degenerate active-shell path tangent")
    result = peak.copy()
    result[active] += 0.01 * tangent / length
    return result


def run_ts(
    campaign: Path, branch_root: Path, rows: list[dict[str, Any]]
) -> dict[str, Any]:
    """Conditionally search for an unrestrained relay TS candidate.

    Args:
        campaign: Campaign with the unchanged w033 native model in ``model``.
        branch_root: Existing relay branch directory.
        rows: Ordered accepted profile rows, including baseline and release.

    Returns:
        Decision or force-accepted candidate awaiting curvature verification.

    Raises:
        ValueError: If scientific inputs, force, or frozen exterior fail.
        FileExistsError: If this search has already been attempted.
    """
    selected = _peak(rows)
    out = branch_root / "03_ts_search"
    journal("intent", f"Create conditional relay TS decision in {out}")
    out.mkdir(exist_ok=False)
    if selected is None:
        decision = {
            "status": "not_launched",
            "reason": "no_force_uphill_interior_maximum",
        }
        write_json(out / "decision.json", decision)
        journal(
            "outcome",
            f"No interior physical-energy maximum; no TS launch in {out}",
            (out / "decision.json",),
        )
        return decision
    before, peak, after = selected
    model = campaign / "model"
    reference_path = model / "inputs/optimized.ncrst"
    reference = NativeState(model, reference_path)
    states = [NativeState(model, Path(row["restart"])) for row in selected]
    if any(not np.array_equal(state.active, reference.active) for state in states):
        raise ValueError("Active atom masks differ")
    mode_xyz = _mode(
        states[0].xyz, states[1].xyz, states[2].xyz, reference.active, reference.xyz
    )
    peak_copy, reference_copy, mode_path = (
        out / name for name in ("peak.ncrst", "reference.ncrst", "mode.ncrst")
    )
    shutil.copy2(Path(peak["restart"]), peak_copy)
    shutil.copy2(reference_path, reference_copy)
    write_restart(mode_path, mode_xyz, states[1].restart)
    mdin = out / "dimer.mdin"
    mdin.write_text(
        _dimer_input(
            reference.mdin_source, reference_copy.resolve(), mode_path.resolve()
        )
    )
    decision = {
        "status": "launching",
        "peak_label": peak["label"],
        "neighbor_labels": [before["label"], after["label"]],
        "peak_physical_kcal_mol": peak["physical_kcal_mol"],
    }
    write_json(out / "decision.json", decision)
    journal(
        "intent",
        f"Run unrestrained dimer candidate at {out}",
        (out / "decision.json", mdin),
    )
    from .nterminal_relay_native import run_native

    result = run_native(
        campaign, out / "native_dimer", peak_copy, mdin.read_text(), restraints=None
    )
    if result.get("accepted") is not True or result.get("native_converged") is not True:
        raise ValueError("Native dimer did not converge")
    candidate = NativeState(model, Path(result["restart"]))
    if not np.array_equal(candidate.active, reference.active) or not np.array_equal(
        candidate.xyz[~reference.active], reference.xyz[~reference.active]
    ):
        raise ValueError("TS candidate changed the frozen exterior")
    energy, gradient = candidate.gradient(
        candidate.xyz, out / "independent_gradient", precision=True
    )
    gradient = np.asarray(gradient)
    if (
        gradient.shape != (int(reference.active.sum()), 3)
        or not np.isfinite(gradient).all()
        or not math.isfinite(energy)
    ):
        raise ValueError("Invalid independent native gradient")
    maximum = float(np.max(np.abs(gradient)) * 0.529177)
    if maximum > 1e-4:
        raise ValueError("Independent TS maximum gradient exceeds 1e-4 Eh/bohr")
    decision.update(
        status="awaiting_curvature_verification",
        native_result=result,
        independent_energy_hartree=energy,
        maximum_gradient_Eh_per_bohr=maximum,
        interpretation="Force-accepted candidate; Hessian index and connectivity are not certified.",
    )
    write_json(out / "decision.json", decision)
    journal(
        "outcome",
        f"TS candidate force accepted; curvature pending at {out}",
        (out / "decision.json",),
    )
    return decision
