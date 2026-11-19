"""Run one explicit capped native DL-FIND L-BFGS relay stage."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7

from .dlfind_lbfgs import DLFindLBFGS
from .native_sander_cli import (
    ANGSTROM_PER_BOHR,
    HARTREE_TO_KCAL_PER_MOL,
    NativeSanderEvaluator,
)
from .relay_stage_units import capped_step, runtime_environment
from .reaction_bias import DistanceCV, HarmonicTarget, ReactionBiasCallback
from .relay_case import (
    CheckedPhysical,
    RelayCase,
    load_case,
    state_record,
    write_json,
    snapshot_stage_inputs,
)
from scripts.prepare_dimer_continuation import write_restart

LOGGER = logging.getLogger(__name__)
LBF_SOURCE = Path(
    "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/libdlfind/src/dlf_lbfgs.F90"
)


def run_stage(
    case_path: Path,
    start_restart: Path,
    output_dir: Path,
    targets: list[HarmonicTarget],
    *,
    max_steps: int,
) -> dict[str, Any]:
    """Run a capped stage with separate process and chemical records.

    Args:
        case_path: Explicit scientific case JSON.
        start_restart: Full-system starting restart.
        output_dir: New stage output directory.
        targets: Analytic terms; empty means zero physical bias.
        max_steps: Nonnegative accepted-step budget.

    Returns:
        Endpath, forceconv, iteration count and evidence location.
    """
    if type(max_steps) is not int or max_steps < 0:
        raise ValueError("max_steps must be a nonnegative integer")
    case_path, start_restart, output_dir = (
        p.resolve() for p in (case_path, start_restart, output_dir)
    )
    model = load_case(case_path)
    restart = Rst7.open(str(start_restart))
    xyz = np.asarray(restart.coordinates, dtype=float).reshape(-1, 3).copy()
    if xyz.shape != model.coordinates.shape or not np.isfinite(xyz).all():
        raise ValueError("Stage starting coordinates mismatch source atom count")
    if not np.array_equal(restart.box, model.restart.box):
        raise ValueError("Stage periodic box differs from immutable case source")
    if not np.array_equal(xyz[~model.active], model.coordinates[~model.active]):
        raise ValueError("Stage frozen exterior differs from immutable case source")
    output_dir.mkdir(parents=True, exist_ok=False)
    old_environment = {
        key: os.environ.get(key)
        for key in (
            "AMBERHOME",
            "LD_LIBRARY_PATH",
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
        )
    }
    amber = Path(model.config["amber_home"])
    environment = runtime_environment(amber)
    os.environ.update(environment)
    try:
        return _run(
            model,
            case_path,
            start_restart,
            output_dir,
            xyz,
            targets,
            max_steps,
            environment,
        )
    except BaseException as error:
        write_json(
            output_dir / "process.json",
            {
                "completed": False,
                "error_type": type(error).__name__,
                "error": str(error),
            },
        )
        write_json(
            output_dir / "scientific_acceptance.json",
            {
                "status": "not_assessed_due_to_process_failure",
                "chemical_acceptance": None,
            },
        )
        raise
    finally:
        for key, value in old_environment.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _run(
    model: RelayCase,
    case_path: Path,
    start: Path,
    output: Path,
    xyz: np.ndarray,
    targets: list[HarmonicTarget],
    max_steps: int,
    environment: dict[str, str],
) -> dict[str, Any]:
    cfg = model.config
    snapshot_stage_inputs(case_path, start, output, model, LBF_SOURCE)
    (output / "one_cycle.mdin").write_text(cfg["mdin"])
    write_json(
        output / "stage_input.json",
        {
            "targets": [t.to_metadata() for t in targets],
            "max_steps": max_steps,
            "environment": environment,
            "memory": cfg["memory"],
            "max_atom_step_A": cfg["max_atom_step_A"],
            "force_gate_kcal_mol_A": cfg["force_gate"],
            "units": {
                "optimizer_coordinates": "bohr",
                "optimizer_gradient": "Eh/bohr",
                "records_gradient": "kcal/mol/A",
            },
        },
    )
    evaluations = output / "evaluations"
    evaluations.mkdir()
    native = NativeSanderEvaluator(
        Path(cfg["topology"]),
        model.restart,
        model.active,
        cfg["mdin"],
        Path(cfg["amber_home"]) / "bin/sander",
        evaluations,
        model.scale,
    )
    if native.counter != 0:
        raise ValueError("New stage native evaluation counter must start at zero")
    effective_targets = targets or [
        HarmonicTarget(
            DistanceCV(model.atom_map["thr_og"] - 1, model.atom_map["carbon"] - 1),
            float(
                np.linalg.norm(
                    xyz[model.atom_map["thr_og"] - 1]
                    - xyz[model.atom_map["carbon"] - 1]
                )
            ),
            0.0,
        )
    ]
    write_json(
        output / "callback_targets.json",
        {
            "unbiased": not bool(targets),
            "zero_spring_adapter": not bool(targets),
            "effective_targets": [t.to_metadata() for t in effective_targets],
        },
    )
    callback = ReactionBiasCallback(
        CheckedPhysical(native, model.active),
        xyz,
        np.flatnonzero(model.active),
        effective_targets,
    )
    x = xyz[model.active].ravel().copy()
    previous = None
    forceconv = False
    with DLFindLBFGS(
        len(x),
        Path(cfg["amber_home"]) / "lib/libdlfind.so",
        LBF_SOURCE,
        memory=cfg["memory"],
    ) as optimizer:
        write_json(output / "optimizer_provenance.json", optimizer.provenance)
        with (output / "accepted_iterations.jsonl").open("x") as records:
            for iteration in range(max_steps + 1):
                energy, gradient = callback(x)
                evaluation = callback.last_evaluation
                assert evaluation is not None
                xyz = callback.coordinates(x)
                if not np.array_equal(
                    xyz[~model.active], model.coordinates[~model.active]
                ):
                    raise ValueError("Frozen exterior changed during stage")
                force = float(
                    np.linalg.norm(gradient.reshape(-1, 3), axis=1).max(initial=0)
                )
                physical_force = float(
                    np.linalg.norm(
                        evaluation.gradient[model.active]
                        - evaluation.bias.gradient[model.active],
                        axis=1,
                    ).max(initial=0)
                )
                max_component = float(np.abs(gradient).max(initial=0))
                forceconv = max_component <= cfg["force_gate"]
                record = {
                    "iteration": iteration,
                    "native_evaluation": native.counter,
                    "physical_energy_kcal_mol": evaluation.physical_energy,
                    "bias_energy_kcal_mol": evaluation.bias.energy,
                    "total_energy_kcal_mol": energy,
                    "max_gradient_kcal_mol_A": max_component,
                    "max_atom_gradient_norm_kcal_mol_A": force,
                    "physical_max_gradient_kcal_mol_A": physical_force,
                    "force_converged": forceconv,
                    "target_values": list(evaluation.bias.values),
                    "state": state_record(model, xyz),
                }
                records.write(json.dumps(record, allow_nan=False) + "\n")
                records.flush()
                if (
                    iteration == 0
                    or iteration % 25 == 0
                    or forceconv
                    or iteration == max_steps
                ):
                    write_restart(
                        output / f"accepted_{iteration:05d}.ncrst", xyz, model.restart
                    )
                if forceconv or iteration == max_steps:
                    break
                xb = x / ANGSTROM_PER_BOHR
                gb = gradient * ANGSTROM_PER_BOHR / HARTREE_TO_KCAL_PER_MOL
                if previous is not None:
                    curvature = float(np.dot(xb - previous[0], gb - previous[1]))
                    if curvature <= 0:
                        LOGGER.info(
                            "Reset native LBFGS history at iteration %d: curvature %.9g",
                            iteration,
                            curvature,
                        )
                        optimizer.restart()
                        with (output / "history_resets.jsonl").open("a") as resets:
                            resets.write(
                                json.dumps(
                                    {
                                        "iteration": iteration,
                                        "secant_curvature_Eh": curvature,
                                        "reason": "nonpositive secant curvature",
                                    }
                                )
                                + "\n"
                            )
                step = optimizer.step(xb, gb)
                if (
                    not np.isfinite(step).all()
                    or np.dot(step, gb) >= 0
                    or not np.any(step)
                ):
                    raise ValueError(
                        "Native LBFGS did not provide a finite nonzero descent direction"
                    )
                previous = (xb.copy(), gb.copy())
                x += capped_step(step, cfg["max_atom_step_A"])
    endpoint = output / "optimized.ncrst"
    write_restart(endpoint, xyz, model.restart)
    result = {
        "endpath": str(endpoint),
        "forceconv": forceconv,
        "iterations": iteration,
        "native_evaluations": native.counter,
        "output_dir": str(output),
    }
    write_json(output / "process.json", {"completed": True, **result})
    write_json(
        output / "scientific_acceptance.json",
        {
            "status": "controller_chemical_review_required",
            "force_converged": forceconv,
            "chemical_acceptance": None,
            "final_state": state_record(model, xyz),
        },
    )
    return result
