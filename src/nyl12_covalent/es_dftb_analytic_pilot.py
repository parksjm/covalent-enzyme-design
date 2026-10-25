"""Explicit decision gate for a bounded pilot with unvalidated numerical forces."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_analytic_pilot(
    engine: dict, force_path: Path, decision_path: Path, inventory_source: Path
) -> tuple[dict, list[Path], dict]:
    """Validate an explicit user-steered analytic-gradient pilot decision.

    Args:
        engine: Source-pinned declared engine configuration.
        force_path: Retained failed numerical force-check report.
        decision_path: Explicit bounded-pilot authorization manifest.
        inventory_source: Independently audited approved QM region config.

    Returns:
        Verified analytic baseline, evidence paths and the retained decision.

    Raises:
        ValueError: Decision scope, identity, physics or finite baseline fails.
    """
    decision = json.loads(decision_path.read_text())
    expected = {
        "mode": "bounded_analytic_gradient_pilot",
        "numerical_force_check_required": False,
        "max_cycles": 25,
        "mm_restraint_weight": 10.0,
        "reaction_coordinate_restraints": False,
        "scope": "Exploratory ES relaxation only; no chemical-model validation",
    }
    if any(
        decision[key] != value
        or (isinstance(value, bool) and decision[key] is not value)
        for key, value in expected.items()
    ):
        raise ValueError("Explicit bounded analytic-pilot scope is invalid")
    if any(
        not isinstance(decision[key], str) or not decision[key].strip()
        for key in ["rationale", "user_steering_reference"]
    ):
        raise ValueError("Explicit rationale and user steering reference required")
    approved = decision["approved_inventory_source"]
    if Path(approved["path"]).resolve() != inventory_source.resolve() or approved[
        "sha256"
    ] != _sha(inventory_source):
        raise ValueError("Pilot decision inventory identity differs from approved case")
    failed = decision["retained_failed_force_check"]
    if Path(failed["path"]).resolve() != force_path.resolve() or failed[
        "sha256"
    ] != _sha(force_path):
        raise ValueError("Retained failed numerical proof identity mismatch")
    force_check = json.loads(force_path.read_text())
    if force_check["accepted"] is not False:
        raise ValueError(
            "Analytic pilot must retain the failed FD result without reinterpretation"
        )
    directory = Path(decision["baseline_directory"])
    result_path, forces_path, request_path = [
        directory / name
        for name in ["result.json", "analytic_forces.npy", "request.json"]
    ]
    if (
        _sha(result_path) != decision["baseline_result_sha256"]
        or _sha(forces_path) != decision["baseline_forces_sha256"]
    ):
        raise ValueError("Analytic baseline artifact identity mismatch")
    baseline, request = (
        json.loads(result_path.read_text()),
        json.loads(request_path.read_text()),
    )
    if request["perturbation"] is not None or baseline["perturbation"] is not None:
        raise ValueError("An unchanged-coordinate analytic baseline is required")
    for key in [
        "topology",
        "coordinates",
        "mm_options",
        "qm_options",
        "expected_inventory",
        "parameter_files",
        "approved_inventory_source",
    ]:
        if request["config"][key] != engine[key]:
            raise ValueError("Analytic baseline source or declared physics mismatch")
    forces = np.load(forces_path, allow_pickle=False)
    box = np.asarray(baseline["box"], dtype=float)
    if (
        not math.isfinite(baseline["energy_kcal_per_mol"])
        or forces.ndim != 2
        or forces.shape[1] != 3
        or not np.isfinite(forces).all()
        or box.shape != (6,)
        or not np.isfinite(box).all()
        or not np.isfinite(baseline["forces"]).all()
    ):
        raise ValueError(
            "Finite full analytic forces, energy and fixed box are required"
        )
    return baseline, [decision_path, result_path, forces_path, request_path], decision
