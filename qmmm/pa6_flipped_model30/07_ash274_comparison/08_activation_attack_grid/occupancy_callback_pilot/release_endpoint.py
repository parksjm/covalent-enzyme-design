"""Test the closest attack endpoint after removing both driving biases."""

from __future__ import annotations

from dataclasses import replace
import json
import logging
from pathlib import Path
import shutil
from typing import Any

import numpy as np
from parmed.amber import Rst7

from nyl12_covalent.native_sander_cli import source_sha256
from nyl12_covalent.occupancy_bias import ActiveBiasCallback
from resume_pilot import ContinuedOptimizer
from runtime_support import (
    ROOT,
    configuration,
    environment,
    load_model,
    save,
    snapshot_sources,
    state_record,
)

LOGGER = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
PILOT = HERE / "pilot"
POINT = PILOT / "point_D"


def classify(xyz: np.ndarray, inventory: dict[str, Any]) -> dict[str, Any]:
    """Describe geometric chemistry without asserting Hessian-verified stability."""
    bound: dict[int, list[int]] = {}
    for proton in inventory["proton_owners"]:
        if proton["nearest_distance_A"] <= 1.25:
            owner = proton["nearest_qm_heavy_1based"]
            bound.setdefault(owner, []).append(proton["hydrogen_1based"])

    def distance(first: int, second: int) -> float:
        return float(np.linalg.norm(xyz[first - 1] - xyz[second - 1]))

    attack = distance(13916, 19920)
    carbonyl = distance(19920, 19921)
    leaving = distance(19920, 19933)
    n_h, o_h = len(bound.get(13905, [])), len(bound.get(13916, []))
    expected_spectators = (
        len(bound.get(20317, [])) == 2
        and len(bound.get(14578, []) + bound.get(14579, [])) == 1
        and not bound.get(14559, []) + bound.get(14560, [])
        and len(bound.get(19933, [])) == 1
    )
    if attack > 2.5 and (n_h, o_h) == (2, 1):
        state = "reactant_E-S_geometry"
    elif attack > 2.5 and (n_h, o_h) == (3, 0):
        state = "proton_transferred_non_covalent_geometry"
    elif (
        1.30 <= attack <= 1.75
        and 1.25 <= carbonyl <= 1.55
        and 1.25 <= leaving <= 1.65
        and (n_h, o_h) == (3, 0)
        and not bound.get(19921)
        and expected_spectators
    ):
        state = "TI1_like_geometry_candidate"
    elif 1.25 <= attack <= 1.65 and leaving > 2.0 and 1.15 <= carbonyl <= 1.30:
        state = "acylenzyme_like_geometry_candidate"
    else:
        state = "other_or_partially_reacted_geometry"
    return {
        "state": state,
        "OG1_C_A": attack,
        "C_O_A": carbonyl,
        "C_N_leaving_A": leaving,
        "N_covalent_H": bound.get(13905, []),
        "OG1_covalent_H": bound.get(13916, []),
        "water_covalent_H": bound.get(20317, []),
        "carbonyl_O_covalent_H": bound.get(19921, []),
        "Asp272_O_covalent_H": bound.get(14559, []) + bound.get(14560, []),
        "Asp274_O_covalent_H": bound.get(14578, []) + bound.get(14579, []),
        "expected_spectator_proton_inventory": expected_spectators,
        "geometric_H_threshold_A": 1.25,
        "classification_is_geometric_not_electronic": True,
        "hessian_minimum_verified": False,
    }


class ReleaseOptimizer(ContinuedOptimizer):
    """Use the same native optimizer and enforce exactly zero added bias."""

    def record(self, y: np.ndarray, label: str) -> dict[str, Any]:
        """Retain inherited accepted-step evidence and verify both springs absent."""
        row = super().record(y, label)
        if row["bias_energy_kcal_mol"] != 0.0:
            raise ValueError("Unbiased release acquired a nonzero bias")
        if row["total_energy_kcal_mol"] != row["physical_energy_kcal_mol"]:
            raise ValueError("Unbiased physical and optimizer energies differ")
        return row


def main() -> None:
    """Run one force-converged candidate's unbiased local optimization in tmux."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    model = load_model()
    settings = environment()
    accepted = json.loads((POINT / "scientific_acceptance.json").read_text())
    if not accepted["biased_force_converged"]:
        raise ValueError("Source D did not meet the original force criterion")
    sources = list(POINT.glob("accepted_*_continuation_final.ncrst"))
    if len(sources) != 1:
        raise ValueError("Expected exactly one final accepted D restart")
    source = sources[0]
    xyz = np.asarray(Rst7.open(str(source)).coordinates).reshape(-1, 3).copy()
    if not np.array_equal(xyz[~model.active], model.coordinates[~model.active]):
        raise ValueError("Source D frozen exterior changed")
    directory = POINT / "unbiased_release"
    directory.mkdir(exist_ok=False)
    paths = model.source_paths + [
        Path(__file__),
        HERE / "resume_pilot.py",
        HERE / "run_pilot.py",
        HERE / "runtime_support.py",
        ROOT / "src/nyl12_covalent/dlfind_lbfgs.py",
        ROOT / "src/nyl12_covalent/occupancy_bias.py",
        ROOT / "src/nyl12_covalent/proton_coordination.py",
        ROOT / "src/nyl12_covalent/native_sander_cli.py",
        ROOT / "scripts/prepare_dimer_continuation.py",
        source,
        POINT / "scientific_acceptance.json",
        PILOT / "continuation_001/protocol.json",
    ]
    hashes = snapshot_sources(directory / "source_snapshot", paths)
    shutil.copy2(source, directory / "starting.ncrst")
    save(
        directory / "protocol.json",
        {
            "hypothesis": "Closest force-converged attack geometry may or may not persist after both driving biases are removed",
            "source_restart": str(source),
            "source_sha256": hashes,
            "source_geometry": classify(xyz, accepted["endpoint"]),
            "bias_kq": 0.0,
            "bias_kr": 0.0,
            "native_ntr": 0,
            "native_nmropt": 0,
            "force_gate_kcal_mol_A": 0.12,
            "max_steps": 1500,
            "movable_atoms": 4662,
            "real_QM_atoms": 102,
            "qm_charge": -1,
            "frozen_exterior": "Unchanged original 12 A complete-residue active selection",
            "method_Hamiltonian_electronic_settings": "Unchanged periodic DFTB3/MM, 0 K electronic temperature",
            "optimizer": "Installed DL-FIND gradient L-BFGS, 0.05 A per-atom step cap",
            "environment": settings,
            "chemical_interpretation": "Source OG1-C about1.854A is not accepted as TI1; release is diagnostic",
            "inherited_output_note": "scientific_acceptance.biased_force_converged field is inherited; with zero springs it measures unbiased force convergence",
            "no_automatic_followup_calculations": True,
        },
    )
    optimizer = ReleaseOptimizer(model, xyz, directory, 2.10, 1.60)
    zero = replace(configuration(model, 2.10, 1.60), k_q=0.0, k_r=0.0)
    optimizer.callback = ActiveBiasCallback(
        optimizer.physical, xyz, optimizer.indices, zero
    )
    (directory / "evaluate.mdin").write_text(model.mdin)
    save(directory / "starting_assessment.json", classify(xyz, accepted["endpoint"]))
    try:
        endpoint, converged = optimizer.run()
    except Exception as error:
        save(
            directory / "process_failure.json",
            {
                "completed": False,
                "exception_type": type(error).__name__,
                "message": str(error),
                "native_evaluations": optimizer.native.counter,
            },
        )
        raise
    final = list(directory.glob("accepted_*_continuation_final.ncrst"))
    if len(final) != 1:
        raise ValueError("Expected exactly one final release restart")
    shutil.copy2(final[0], directory / "optimized.ncrst")
    inventory = state_record(model, endpoint)
    assessment = classify(endpoint, inventory)
    assessment.update(
        {
            "unbiased_force_converged": converged,
            "both_driving_biases_removed": True,
            "endpoint_inventory": inventory,
            "source_and_final_sha256": {
                str(p): source_sha256(p) for p in [source, final[0]]
            },
            "full_stationary_point_verification_pending": True,
        }
    )
    save(directory / "release_assessment.json", assessment)
    LOGGER.info(
        "Unbiased endpoint: %s; force_converged=%s", assessment["state"], converged
    )


if __name__ == "__main__":
    main()
