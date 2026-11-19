"""Prove energy/force coupling through the actual optimizer, not a native CV hook."""

from __future__ import annotations

from dataclasses import replace
import json
import logging
from pathlib import Path
import re
import sys

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

import numpy as np  # noqa: E402
from parmed.amber import Rst7  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

from nyl12_covalent.native_sander_cli import NativeSanderEvaluator  # noqa: E402
from nyl12_covalent.occupancy_bias import ActiveBiasCallback  # noqa: E402
from nyl12_covalent.proton_coordination import coordination_transfer  # noqa: E402
from runtime_support import (  # noqa: E402
    AMBER,
    SOURCE,
    CheckedPhysical,
    configuration,
    environment,
    load_model,
    save,
    snapshot_sources,
    state_record,
)

LOGGER = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent


def main() -> None:
    """Run exactly two optimizer iterations with a known nonzero added bias."""
    logging.basicConfig(level=logging.INFO)
    model = load_model()
    folder = HERE / "callback_proof"
    folder.mkdir(exist_ok=False)
    native_dir = folder / "native"
    native_dir.mkdir()
    restart = Rst7.open(str(SOURCE / "optimized.ncrst"))
    xyz = np.asarray(restart.coordinates, dtype=float).reshape(-1, 3).copy()
    transfer = coordination_transfer(xyz, 13915, 13904, model.hydrogens)
    distance = float(np.linalg.norm(xyz[13915] - xyz[19919]))
    config = replace(
        configuration(model, transfer.value + 0.2, distance), k_q=10.0, k_r=0.0
    )
    settings = environment()
    sources = snapshot_sources(
        folder / "source_snapshot",
        model.source_paths
        + [
            Path(__file__),
            HERE / "runtime_support.py",
            ROOT / "src/nyl12_covalent/occupancy_bias.py",
            ROOT / "src/nyl12_covalent/proton_coordination.py",
            ROOT / "src/nyl12_covalent/native_sander_cli.py",
        ],
    )
    save(
        folder / "protocol.json",
        {
            "scope": "two-iteration functional proof, not stationary-point optimization",
            "source_sha256": sources,
            "environment": settings,
            "source_restart": str(SOURCE / "optimized.ncrst"),
            "expected_initial_bias_kcal_mol": 0.2,
            "q_initial": transfer.value,
            "q_target": config.q_target,
            "k_q": 10.0,
            "active_atoms": int(model.active.sum()),
            "coordinate_scale_A": 0.1,
        },
    )
    native = NativeSanderEvaluator(
        model.topology_path,
        restart,
        model.active,
        model.mdin,
        AMBER / "bin/sander",
        native_dir,
        model.scale,
    )
    checked = CheckedPhysical(native, model.active)
    physical_result: dict = {}

    def physical(coords: np.ndarray) -> tuple[float, np.ndarray]:
        energy, gradient = checked(coords)
        physical_result.update(energy=energy, gradient=gradient.copy())
        return energy, gradient

    active = np.flatnonzero(model.active)
    callback = ActiveBiasCallback(physical, xyz, active, config)
    records = []
    origin = []

    def objective(y: np.ndarray) -> tuple[float, np.ndarray]:
        total, gradient = callback(y * 0.1)
        value = callback.last_evaluation
        if value is None:
            raise ValueError("Optimizer callback omitted bias assessment")
        expected_gradient = value.bias.gradient
        error = np.max(
            np.abs(value.gradient - physical_result["gradient"] - expected_gradient)
        )
        rebuilt = callback.coordinates(y * 0.1)
        if not np.array_equal(rebuilt[~model.active], xyz[~model.active]):
            raise ValueError("Frozen outer coordinates changed")
        records.append(
            {
                "evaluation": len(records) + 1,
                "physical_kcal_mol": value.physical_energy,
                "total_kcal_mol": total,
                "bias_kcal_mol": value.bias.energy,
                "q": value.bias.q,
                "added_gradient_max_kcal_mol_A": float(
                    np.max(np.abs(expected_gradient))
                ),
                "gradient_composition_max_error": float(error),
                "native_evidence": native.evidence[-1],
            }
        )
        with (folder / "evaluations.jsonl").open("a") as stream:
            stream.write(json.dumps(records[-1]) + "\n")
        if not origin:
            origin.append(total)
        return total - origin[0], gradient * 0.1

    accepted = []

    def iteration(y: np.ndarray) -> None:
        accepted.append(state_record(model, callback.coordinates(y * 0.1)))

    result = minimize(
        objective,
        xyz[active].ravel() / 0.1,
        method="L-BFGS-B",
        jac=True,
        callback=iteration,
        options={"maxiter": 2, "maxfun": 20, "maxls": 10, "ftol": 0.0, "gtol": 0.001},
    )
    first = records[0]
    native_energies = re.findall(
        r"EPtot\s*=\s*([-+\d.EeDd]+)", (native_dir / "eval_00001.mdout").read_text()
    )
    physical_energy_error = abs(
        first["physical_kcal_mol"] - float(native_energies[0].replace("D", "E"))
    )
    gates = {
        "known_bias_applied": abs(first["bias_kcal_mol"] - 0.2) < 1e-10,
        "energy_composition": abs(
            first["total_kcal_mol"] - first["physical_kcal_mol"] - 0.2
        )
        < 1e-8,
        "nonzero_added_force": first["added_gradient_max_kcal_mol_A"] > 1e-6,
        "analytic_force_composition": all(
            r["gradient_composition_max_error"] < 1e-9 for r in records
        ),
        "native_energy_unit_match": physical_energy_error < 0.01,
        "optimizer_exercised": len(records) > 1 and len(accepted) > 0,
        "native_input_geometry_match": checked.initial_geometry_error_A is not None
        and checked.initial_geometry_error_A <= 2e-6,
    }
    report = {
        "accepted": all(gates.values()),
        "gates": gates,
        "native_energy_conversion_error_kcal_mol": physical_energy_error,
        "optimizer_message": str(result.message),
        "optimizer_iterations": int(result.nit),
        "native_evaluations": len(records),
        "initial_bias_kcal_mol": first["bias_kcal_mol"],
        "accepted_optimizer_states": accepted,
        "stationary_point_claim": False,
    }
    save(folder / "scientific_acceptance.json", report)
    if not report["accepted"]:
        raise RuntimeError(f"Biased optimizer proof failed: {gates}")
    LOGGER.info(
        "Optimizer proof passed: %.6f kcal/mol added; %d native evaluations",
        first["bias_kcal_mol"],
        len(records),
    )


if __name__ == "__main__":
    main()
