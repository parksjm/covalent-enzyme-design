"""Continue the pilot in place using installed DL-FIND gradient-based L-BFGS."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7

from nyl12_covalent.dlfind_lbfgs import DLFindLBFGS
from run_pilot import FORCE_GATE, SCALE_A, StageOptimizer
from runtime_support import (
    AMBER,
    POINTS,
    ROOT,
    environment,
    load_model,
    save,
    snapshot_sources,
    state_record,
)

LOGGER = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
PILOT = HERE / "pilot"
SOURCE = ROOT / (
    "../code/ambertools26/build_assets/src/ambertools26_src/"
    "AmberTools/src/libdlfind/src/dlf_lbfgs.F90"
)
BOHR_A = 0.529177
HARTREE_KCAL = 627.5094740631
MAX_STEPS = 1500
MAX_ATOM_STEP_A = 0.05


def archive_terminal_records() -> None:
    """Preserve failed terminal records before reopening the same logical pilot."""
    pairs = [
        (PILOT / "completion.json", PILOT / "history/line_search_stop/completion.json"),
        (
            PILOT / "point_A/process.json",
            PILOT / "point_A/history/line_search_stop/process.json",
        ),
        (
            PILOT / "point_A/scientific_acceptance.json",
            PILOT / "point_A/history/line_search_stop/scientific_acceptance.json",
        ),
    ]
    for source, destination in pairs:
        if not source.is_file() or destination.exists():
            raise ValueError(f"Unexpected archive state for {source}")
    for source, destination in pairs:
        destination.parent.mkdir(parents=True, exist_ok=True)
        source.rename(destination)


class ContinuedOptimizer(StageOptimizer):
    """Retain original bias, evaluator and recorder; replace step selection only."""

    def run(self) -> tuple[np.ndarray, bool]:
        """Take capped native L-BFGS steps without an energy line search."""
        y = self.coordinates_initial[self.indices].ravel() / SCALE_A
        self.record(y, "continuation_initial")
        previous_x: np.ndarray | None = None
        previous_g: np.ndarray | None = None
        events = self.directory / "continuation_optimizer_events.jsonl"
        with DLFindLBFGS(
            n_variables=y.size,
            library_path=AMBER / "lib/libdlfind.so",
            source_path=SOURCE.resolve(strict=True),
            memory=7,
        ) as optimizer:
            provenance: dict[str, Any] = optimizer.provenance
            for _ in range(MAX_STEPS):
                _, gy = self.objective(y)
                gradient = gy / SCALE_A
                if float(np.max(np.abs(gradient))) <= FORCE_GATE:
                    break
                x = y * SCALE_A / BOHR_A
                g = gradient * BOHR_A / HARTREE_KCAL
                if previous_x is not None and previous_g is not None:
                    secant = float(np.dot(x - previous_x, g - previous_g))
                    if secant <= 0 or not np.any(g != previous_g):
                        with events.open("a") as stream:
                            stream.write(
                                json.dumps(
                                    {
                                        "accepted_iteration": self.accepted,
                                        "event": "explicit_nonpositive_curvature_history_restart",
                                        "secant_bohr_hartree": secant,
                                    }
                                )
                                + "\n"
                            )
                        optimizer.restart()
                        LOGGER.info(
                            "Explicit L-BFGS history restart: secant=%g", secant
                        )
                step = optimizer.step(x, g)
                if float(np.dot(step, g)) >= 0:
                    raise ValueError(
                        "Native L-BFGS did not propose a descent direction"
                    )
                step_A = step.reshape(-1, 3) * BOHR_A
                largest = float(np.max(np.linalg.norm(step_A, axis=1)))
                if not np.isfinite(largest) or largest == 0:
                    raise ValueError("Native L-BFGS produced an invalid step")
                factor = min(1.0, MAX_ATOM_STEP_A / largest)
                previous_x, previous_g = x.copy(), g.copy()
                y = y + factor * step_A.ravel() / SCALE_A
                self.accepted_step(y)
            final = self.record(y, "continuation_final")
        converged = final["max_active_gradient_kcal_mol_A"] <= FORCE_GATE
        coordinates = self.callback.coordinates(y * SCALE_A)
        save(
            self.directory / "process.json",
            {
                "completed": True,
                "native_process_and_SCC_passed": True,
                "optimizer": "Installed DL-FIND gradient-based L-BFGS step routine",
                "optimizer_provenance": provenance,
                "energy_line_search": False,
                "max_atom_step_A": MAX_ATOM_STEP_A,
                "accepted_iteration_final": self.accepted,
                "native_evaluations_final": self.native.counter,
                "termination": "force_converged" if converged else "iteration_budget",
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
    """Validate the failed A checkpoint, preserve it, and continue A–D serially."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    model = load_model()
    settings = environment()
    completion = json.loads((PILOT / "completion.json").read_text())
    process = json.loads((PILOT / "point_A/process.json").read_text())
    if completion != {
        "status": "incomplete",
        "stopped_after": "A",
        "automatic_extension": False,
    }:
        raise ValueError("Unexpected prior pilot completion state")
    if process["optimizer_status"] != 2 or process["accepted_iterations"] != 158:
        raise ValueError("Unexpected original failed point A")
    initial = PILOT / "point_A/accepted_00158_final.ncrst"
    coordinates = np.asarray(Rst7.open(str(initial)).coordinates).reshape(-1, 3).copy()
    if not np.array_equal(coordinates[~model.active], model.coordinates[~model.active]):
        raise ValueError("Frozen exterior differs from original source")
    if any((PILOT / f"point_{name}").exists() for name in "BCD"):
        raise ValueError("Continuation expects B–D to be unstarted")
    metadata = PILOT / "continuation_001"
    metadata.mkdir(exist_ok=False)
    sources = model.source_paths + [
        Path(__file__).resolve(),
        HERE / "run_pilot.py",
        HERE / "runtime_support.py",
        ROOT / "src/nyl12_covalent/dlfind_lbfgs.py",
        SOURCE.resolve(strict=True),
        ROOT / "src/nyl12_covalent/occupancy_bias.py",
        ROOT / "src/nyl12_covalent/proton_coordination.py",
        ROOT / "src/nyl12_covalent/native_sander_cli.py",
        ROOT / "scripts/prepare_dimer_continuation.py",
        initial,
        PILOT / "protocol.json",
        PILOT / "completion.json",
        PILOT / "point_A/process.json",
        PILOT / "point_A/scientific_acceptance.json",
        HERE / "line_search_diagnosis/audit.json",
    ]
    hashes = snapshot_sources(metadata / "source_snapshot", sources)
    save(
        metadata / "protocol.json",
        {
            "starting_restart": str(initial),
            "source_sha256": hashes,
            "environment": settings,
            "points": POINTS,
            "change": "Replace energy-Wolfe line search with installed DL-FIND gradient-based L-BFGS",
            "reason": "Native MM10A cutoff jump exceeds final predicted energy decrease",
            "geometry_restraints_QM_charge_shell_Hamiltonian_force_gate": "unchanged",
            "memory": 7,
            "max_steps_each": MAX_STEPS,
            "max_atom_step_A": MAX_ATOM_STEP_A,
            "energy_line_search": False,
            "units": "Native step routine uses bohr and hartree/bohr; callback kcal/mol and angstrom",
            "history_restart": "Explicitly recorded if secant curvature nonpositive",
            "no_new_dependency": True,
        },
    )
    archive_terminal_records()
    for name, q_target, r_target in POINTS:
        directory = PILOT / f"point_{name}"
        if name != "A":
            directory.mkdir()
            save(
                directory / "stage.json",
                {
                    "name": name,
                    "q_target": q_target,
                    "r_target_A": r_target,
                    "start": state_record(model, coordinates),
                    "continuation_protocol": str(metadata / "protocol.json"),
                },
            )
            snapshot_sources(
                directory / "source_snapshot", [metadata / "protocol.json"]
            )
            (directory / "evaluate.mdin").write_text(model.mdin)
        optimizer = ContinuedOptimizer(
            model, coordinates, directory, q_target, r_target
        )
        if name == "A":
            optimizer.accepted = 158
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
                PILOT / "completion.json",
                {
                    "status": "incomplete",
                    "stopped_after": name,
                    "automatic_extension": False,
                    "continuation": "001",
                },
            )
            return
    save(
        PILOT / "completion.json",
        {
            "status": "four_candidates_force_converged",
            "chemical_identity_validated": False,
            "continuation": "001",
        },
    )


if __name__ == "__main__":
    main()
