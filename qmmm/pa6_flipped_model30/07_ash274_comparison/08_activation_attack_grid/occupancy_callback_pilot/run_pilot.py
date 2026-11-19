"""Prepare or run the explicitly bounded serial occupancy-bias diagnostic pilot."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize

from nyl12_covalent.native_sander_cli import NativeSanderEvaluator
from nyl12_covalent.occupancy_bias import ActiveBiasCallback
from scripts.prepare_dimer_continuation import write_restart
from runtime_support import (
    AMBER,
    POINTS,
    ROOT,
    CheckedPhysical,
    Model,
    configuration,
    environment,
    load_model,
    preflight,
    save,
    snapshot_sources,
    state_record,
)

LOGGER = logging.getLogger(__name__)
SCALE_A = 0.1
FORCE_GATE = 0.12
OPTIONS = {
    "gtol": SCALE_A * FORCE_GATE,
    "ftol": 0.0,
    "maxiter": 1500,
    "maxfun": 3000,
    "maxls": 20,
    "maxcor": 10,
}


class StageOptimizer:
    """Scaled optimizer objective and accepted-iteration scientific recorder."""

    def __init__(
        self,
        model: Model,
        coordinates: np.ndarray,
        directory: Path,
        q_target: float,
        r_target: float,
    ) -> None:
        """Create a fresh native evaluator for one bounded pilot point."""
        self.model, self.directory = model, directory
        self.coordinates_initial = coordinates.copy()
        self.indices = np.flatnonzero(model.active)
        self.native = NativeSanderEvaluator(
            model.topology_path,
            model.restart,
            model.active,
            model.mdin,
            AMBER / "bin/sander",
            directory,
            model.scale,
        )
        self.physical = CheckedPhysical(self.native, model.active)
        self.callback = ActiveBiasCallback(
            self.physical,
            coordinates,
            self.indices,
            configuration(model, q_target, r_target),
        )
        self.accepted: int = 0
        self.last_y: np.ndarray | None = None
        self.last_value: tuple[float, np.ndarray] | None = None

    def objective(self, y: np.ndarray) -> tuple[float, np.ndarray]:
        """Return kcal energy and its derivative with respect to x/0.1 Å."""
        if self.last_y is None or not np.array_equal(y, self.last_y):
            energy, gradient = self.callback(y * SCALE_A)
            self.last_y = y.copy()
            self.last_value = energy, gradient * SCALE_A
        if self.last_value is None:
            raise RuntimeError("Objective result unavailable")
        return self.last_value

    def record(self, y: np.ndarray, label: str) -> dict[str, Any]:
        """Record an accepted geometry separately from line-search evaluations."""
        self.objective(y)
        evaluation = self.callback.last_evaluation
        if evaluation is None:
            raise RuntimeError("Missing accepted-geometry evaluation")
        xyz = self.callback.coordinates(y * SCALE_A)
        if not np.array_equal(
            xyz[~self.model.active], self.coordinates_initial[~self.model.active]
        ):
            raise ValueError("Frozen outer coordinates changed")
        bias = evaluation.bias
        record = {
            "accepted_iteration": self.accepted,
            "label": label,
            "native_evaluations": self.native.counter,
            "physical_energy_kcal_mol": evaluation.physical_energy,
            "total_energy_kcal_mol": evaluation.energy,
            "bias_energy_kcal_mol": bias.energy,
            "bias_q_kcal_mol": bias.energy_q,
            "bias_r_kcal_mol": bias.energy_r,
            "CN_N": bias.cn_acceptor,
            "CN_OG1": bias.cn_donor,
            "q": bias.q,
            "r_A": bias.r,
            "q_target_lag": bias.q - self.callback.config.q_target,
            "r_target_lag_A": bias.r - self.callback.config.r_target,
            "max_active_gradient_kcal_mol_A": float(
                np.max(np.abs(evaluation.gradient[self.model.active]))
            ),
            "frozen_outer_exact": True,
        }
        with (self.directory / "accepted_iterations.jsonl").open("a") as stream:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
        if label != "accepted" or self.accepted % 25 == 0:
            name = f"accepted_{self.accepted:05d}_{label}.ncrst"
            write_restart(self.directory / name, xyz, self.model.restart)
        LOGGER.info(
            "%s accepted=%d q=%.6f r=%.6f maxgrad=%.6f",
            self.directory.name,
            self.accepted,
            bias.q,
            bias.r,
            record["max_active_gradient_kcal_mol_A"],
        )
        return record

    def accepted_step(self, y: np.ndarray) -> None:
        """Record optimizer-accepted steps only; trial evaluations remain native I/O."""
        self.accepted += 1
        self.record(y, "accepted")

    def run(self) -> tuple[np.ndarray, bool]:
        """Run once and distinguish optimizer termination from force acceptance."""
        y0 = self.coordinates_initial[self.indices].ravel() / SCALE_A
        self.record(y0, "initial")
        result = minimize(
            self.objective,
            y0,
            method="L-BFGS-B",
            jac=True,
            callback=self.accepted_step,
            options=OPTIONS,
        )
        final = self.record(result.x, "final")
        converged = final["max_active_gradient_kcal_mol_A"] <= FORCE_GATE
        coordinates = self.callback.coordinates(result.x * SCALE_A)
        save(
            self.directory / "process.json",
            {
                "completed": True,
                "native_process_and_SCC_passed": True,
                "optimizer_success": bool(result.success),
                "optimizer_status": int(result.status),
                "optimizer_message": str(result.message),
                "accepted_iterations": int(result.nit),
                "optimizer_function_calls": int(result.nfev),
                "native_evaluations": self.native.counter,
                "initial_native_mapping_max_error_A": self.physical.initial_geometry_error_A,
            },
        )
        save(
            self.directory / "scientific_acceptance.json",
            {
                "biased_force_converged": converged,
                "pilot_stage_accepted": converged,
                "status": "biased_force_converged" if converged else "incomplete",
                "force_gate_kcal_mol_A": FORCE_GATE,
                "max_active_gradient_kcal_mol_A": final[
                    "max_active_gradient_kcal_mol_A"
                ],
                "coordinate_target_lag_is_monitor_only": True,
                "chemical_identity_validated": False,
                "endpoint": state_record(self.model, coordinates),
            },
        )
        return coordinates, converged


def main() -> None:
    """Validate inputs or create one new pilot directory and run serial stages."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    model = load_model()
    report = preflight(model)
    if args.preflight:
        LOGGER.info("%s", json.dumps(report, allow_nan=False))
        return
    if args.output is None:
        parser.error("--output is required for execution")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    settings = environment()
    scripts = [
        Path(__file__).resolve(),
        Path(__file__).with_name("runtime_support.py").resolve(),
        ROOT / "src/nyl12_covalent/occupancy_bias.py",
        ROOT / "src/nyl12_covalent/proton_coordination.py",
        ROOT / "src/nyl12_covalent/native_sander_cli.py",
        ROOT / "scripts/prepare_dimer_continuation.py",
    ]
    hashes = snapshot_sources(output / "source_snapshot", model.source_paths + scripts)
    save(output / "preflight.json", report)
    save(
        output / "protocol.json",
        {
            "purpose": "Four serial constrained-minimization diagnostic candidates; not MFEP/PMF",
            "method": "Native Sander DFTB3/MM, source QM/PME Hamiltonian, electronic temperature 0 K",
            "method_validation": "No DFT-method validation assumed; source import is not certification",
            "protonation": "Original accepted ASP272/ASH274 neutral-Thr water-QM source; all 54 real H mobile",
            "coordinate": "q=CN_N(all54H)-CN_OG1(all54H); independent OG1-substrate C distance",
            "chemical_caution": "qCN alone does not establish chemical identity; inspect all proton owners and original bonds",
            "points": POINTS,
            "springs": {"kq_kcal_mol_CN2": 100.0, "kr_kcal_mol_A2": 200.0},
            "optimizer": {
                "method": "L-BFGS-B",
                "options": OPTIONS,
                "x_scale_A": SCALE_A,
                "extra_bounds": None,
                "extra_position_restraints": None,
            },
            "force_gate_kcal_mol_A": FORCE_GATE,
            "environment": settings,
            "source_sha256": hashes,
            "frozen_atoms": "All atoms outside explicit source 4662-atom active selection",
            "restraints": "Only two stated CV harmonics; no identity hard restraint",
            "continuation": "A starts grid q_004_attack_001; next stage only after force/SCC gates",
        },
    )
    coordinates = model.coordinates.copy()
    for name, q_target, r_target in POINTS:
        directory = output / f"point_{name}"
        directory.mkdir()
        snapshot_sources(directory / "source_snapshot", model.source_paths + scripts)
        (directory / "evaluate.mdin").write_text(model.mdin)
        save(
            directory / "stage.json",
            {
                "name": name,
                "q_target": q_target,
                "r_target_A": r_target,
                "start": state_record(model, coordinates),
            },
        )
        optimizer = StageOptimizer(model, coordinates, directory, q_target, r_target)
        try:
            coordinates, accepted = optimizer.run()
        except Exception as error:
            save(
                directory / "process_failure.json",
                {
                    "completed": False,
                    "exception_type": type(error).__name__,
                    "message": str(error),
                    "native_evaluations_attempted": optimizer.native.counter,
                },
            )
            save(
                directory / "scientific_failure.json",
                {
                    "pilot_stage_accepted": False,
                    "status": "incomplete_process_or_SCC_failure",
                },
            )
            raise
        if not accepted:
            save(
                output / "completion.json",
                {
                    "status": "incomplete",
                    "stopped_after": name,
                    "automatic_extension": False,
                },
            )
            return
    save(
        output / "completion.json",
        {
            "status": "four_candidates_force_converged",
            "chemical_identity_validated": False,
        },
    )


if __name__ == "__main__":
    main()
