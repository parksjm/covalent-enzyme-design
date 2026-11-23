"""Native repeat/parity/step gates and full-shell minimum verification."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.relay_minimum import (
    journal,
    probe_columns,
    proton_inventory,
    verify_frozen,
)
from nyl12_covalent.stationary_hessian import full_hessian, worker_column
from nyl12_covalent.stationary_native import NativeState, write_json


CAVEAT = (
    "Full 12495-DOF active-shell finite-difference inertia at fixed exterior. "
    "The 2*skew threshold is a numerical heuristic, not a rigorous operator-error "
    "bound; selected-column stability does not bound all truncation error. "
    "A pass is numerical local-curvature evidence within these qualifications, "
    "not an analytic proof or chemical/product assignment."
)


def verify_prerequisites(output: Path) -> None:
    """Reject changes to any original method or installed parameter artifact."""
    import json

    from nyl12_covalent.stationary_native import digest

    for name, expected in json.loads(
        (output / "prerequisites.json").read_text()
    ).items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Original prerequisite changed: {name}")


def gradient_metrics(
    energies: list[float], gradients: list[np.ndarray]
) -> dict[str, Any]:
    """Check exact precise repeats, old/new parity and original force tolerance.

    Args:
        energies: Original serialization, precise repeat one and precise repeat two.
        gradients: Corresponding complete active gradients in Eh/Angstrom.

    Returns:
        Independent numerical acceptance metrics.
    """
    if len(energies) != 3 or len(gradients) != 3:
        raise ValueError("Require original serialization and two precise repeats")
    if not np.isfinite(energies).all() or not all(
        np.isfinite(g).all() for g in gradients
    ):
        raise ValueError("Nonfinite gradient/energy evidence")
    if any(g.shape != gradients[0].shape for g in gradients):
        raise ValueError("Gradient dimensions differ")
    exact = energies[1] == energies[2] and np.array_equal(gradients[1], gradients[2])
    parity = float(np.max(np.abs(gradients[0] - gradients[1])))
    bound = 0.50001e-7 / 0.529177**2
    force = float(np.max(np.abs(gradients[1])) * 0.529177)
    difference = abs(energies[0] - energies[1])
    return {
        "passed": exact and parity <= bound and difference <= 1e-9 and force <= 1e-5,
        "precise_repeats_exact": exact,
        "serialization_parity_max_difference": parity,
        "old_rounding_bound_Eh_per_A": bound,
        "energy_difference_Eh": difference,
        "max_gradient_Eh_per_bohr": force,
        "force_tolerance_Eh_per_bohr": 1e-5,
    }


def column_metrics(coarse: np.ndarray, fine: np.ndarray) -> dict[str, Any]:
    """Apply the declared relative step-stability engineering screen."""
    if coarse.shape != fine.shape or not (
        np.isfinite(coarse).all() and np.isfinite(fine).all()
    ):
        raise ValueError("Invalid finite-difference columns")
    norm = float(np.linalg.norm(fine))
    if norm == 0:
        raise ValueError("Zero probe column leaves relative stability undefined")
    difference = float(np.linalg.norm(coarse - fine))
    return {
        "passed": difference / norm <= 1e-4,
        "difference_norm": difference,
        "relative_difference": difference / norm,
        "engineering_relative_gate": 1e-4,
        "limitation": "Selected-column screen only; not a full-matrix error bound",
    }


def precision_preflight(frozen: Path, output: Path, workers: int) -> dict[str, Any]:
    """Block full Hessian until exact repeats/parity and selected steps accept.

    Args:
        frozen: Immutable original Hamiltonian, active mask, restart and library.
        output: Single minimum_verification root.
        workers: At most eight lower-priority single-CPU derivative workers.

    Returns:
        Precision gate retaining every underlying native evaluation.
    """
    verify_frozen(frozen)
    verify_prerequisites(output)
    state = NativeState(frozen)
    if len(state.indices) * 3 != 12495:
        raise ValueError("Full original 12495 Cartesian degrees required")
    stage = output / "02_gradient_repeats"
    stage.mkdir()
    journal(
        "intent",
        "Evaluate original serialization and two exact precise native gradients",
    )
    samples = [
        state.gradient(state.xyz, stage / name, precision=precise)
        for name, precise in (
            ("original", False),
            ("precise_one", True),
            ("precise_two", True),
        )
    ]
    repeats = gradient_metrics([s[0] for s in samples], [s[1] for s in samples])
    repeats["proton_inventory"] = proton_inventory(state)
    write_json(stage / "scientific_acceptance.json", repeats)
    journal(
        "outcome",
        f"Native repeat/parity/force gate passed={repeats['passed']}",
        (stage / "scientific_acceptance.json",),
    )
    if not repeats["passed"]:
        return {"passed": False, "reason": "Native repeat/parity/force gate failed"}
    columns = probe_columns(state)
    stage = output / "03_step_probes"
    stage.mkdir()
    for label in ("h_0p001A", "h_0p0005A"):
        (stage / label).mkdir()
    journal(
        "intent", f"Compute {len(columns)} full-length central columns at h=.001/.0005A"
    )
    tasks = [
        (frozen, stage / label, state.restart_path, column, step)
        for label, step in (("h_0p001A", 0.001), ("h_0p0005A", 0.0005))
        for column in columns
    ]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        list(executor.map(worker_column, tasks))
    results = {}
    for column in columns:
        arrays = [
            np.load(stage / label / "columns" / f"column_{column:05d}" / "column.npy")
            for label in ("h_0p001A", "h_0p0005A")
        ]
        results[str(column)] = column_metrics(*arrays)
    acceptance = {
        "passed": all(r["passed"] for r in results.values()),
        "dimension": 12495,
        "columns": results,
        "step_A": [0.001, 0.0005],
        "meaning": "Bounded engineering screen; no rigorous global error bound",
    }
    write_json(stage / "scientific_acceptance.json", acceptance)
    journal(
        "outcome",
        f"Selected full-length step probes passed={acceptance['passed']}",
        (stage / "scientific_acceptance.json",),
    )
    return acceptance


def verify_minimum(frozen: Path, output: Path, workers: int) -> dict[str, Any]:
    """Classify the complete fixed-exterior active-shell numerical inertia.

    Args:
        frozen: Immutable NativeState root.
        output: Existing single minimum_verification root.
        workers: Declared positive maximum of eight native single-CPU workers.

    Returns:
        Literal full Hessian result with explicit scope and limitations.
    """
    if not 1 <= workers <= 8:
        raise ValueError("Require 1..8 single-CPU workers")
    preflight = precision_preflight(frozen, output, workers)
    if not preflight["passed"]:
        result = {
            "numerically_verified": False,
            "status": "inconclusive_preflight",
            "preflight": preflight,
            "caveat": CAVEAT,
        }
    else:
        verify_frozen(frozen)
        verify_prerequisites(output)
        journal(
            "intent", "Compute complete 12495-DOF Hessian h=.0005A expected_negative=0"
        )
        inertia = full_hessian(
            frozen,
            frozen / "inputs/optimized.ncrst",
            output / "04_full_hessian",
            workers,
            0.0005,
            expected_negative=0,
        )
        result = {
            "numerically_verified": bool(inertia["passed"]),
            "status": "numerical_local_minimum"
            if inertia["passed"]
            else "not_verified",
            "literal_inertia": inertia,
            "caveat": CAVEAT,
        }
        journal(
            "outcome",
            f"Full active-shell Hessian numerical acceptance={inertia['passed']}",
            (output / "04_full_hessian/scientific_acceptance.json",),
        )
    write_json(output / "scientific_acceptance.json", result)
    return result
