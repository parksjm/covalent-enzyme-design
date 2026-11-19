"""Audit retained native calls; never invoke the QM/MM engine."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from itertools import product
import sys
from pathlib import Path
from typing import Any

import numpy as np
from parmed import unit
from parmed.amber import Rst7
from parmed.geometry import box_lengths_and_angles_to_vectors
from scipy.spatial import cKDTree

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = ROOT / (
    "qmmm/pa6_flipped_model30/07_ash274_comparison/"
    "08_activation_attack_grid/occupancy_callback_pilot"
)
sys.path[:0] = [str(ROOT / "src"), str(ROOT), str(BASE)]
from runtime_support import configuration, load_model  # noqa: E402
from nyl12_covalent.occupancy_bias import evaluate_bias  # noqa: E402

OUT = BASE / "line_search_diagnosis"
CALLS = BASE / "pilot/point_A"
EH = 627.5094740631
SCALE = -EH / 0.529177**2


def digest(path: Path) -> str:
    """Return a SHA256 identity for retained evidence."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cutoff_state(xyz: np.ndarray, box: np.ndarray, qm: np.ndarray) -> dict[str, Any]:
    """Describe real-QM/MM 10 Å proximity using 27 triclinic images."""
    cell = np.asarray(
        box_lengths_and_angles_to_vectors(*box).value_in_unit(unit.angstroms),
        dtype=float,
    )
    frac = xyz @ np.linalg.inv(cell)
    wrapped = (frac % 1.0) @ cell
    shifts = np.asarray(list(product((-1, 0, 1), repeat=3))) @ cell
    images = (wrapped[qm][None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    sparse = cKDTree(images).sparse_distance_matrix(
        cKDTree(wrapped), 10.002, output_type="ndarray"
    )
    nearest = np.full(len(xyz), np.inf)
    np.minimum.at(nearest, sparse["j"], sparse["v"])
    nearest[qm] = np.inf
    members = set(np.flatnonzero(nearest <= 10.0).tolist())
    boundary = {
        int(i): float(nearest[i])
        for i in np.flatnonzero(np.abs(nearest - 10.0) <= 0.002)
    }
    return {"members": members, "boundary": boundary}


def main() -> None:
    """Write scalar steps, endpoint comparisons and hashes for the last 60 calls."""
    model = load_model()
    active = np.flatnonzero(model.active)
    config = configuration(model, 2.1, 3.125178400033067)
    records = [
        json.loads(line)
        for line in (CALLS / "native_evidence.jsonl").read_text().splitlines()
    ]
    records = records[-60:]
    hashes = {
        str(path): digest(path)
        for path in (
            Path(__file__),
            BASE / "runtime_support.py",
            ROOT / "src/nyl12_covalent/occupancy_bias.py",
            ROOT / "src/nyl12_covalent/proton_coordination.py",
            ROOT / "src/nyl12_covalent/native_sander_cli.py",
            model.topology_path,
            CALLS / "native_evidence.jsonl",
            CALLS / "eval_00284.mdin",
        )
    }
    states = []
    rows = []
    previous = None
    cutoff_events = []
    for record in records:
        number = record["evaluation"]
        stem = f"eval_{number:05d}"
        restart_path = CALLS / f"{stem}.input.ncrst"
        force_path = CALLS / f"{stem}.path_force.xyz"
        for path in (restart_path, force_path):
            actual = digest(path)
            if actual != record["artifacts_sha256"][path.name]:
                raise ValueError(f"Retained native evidence hash mismatch: {path}")
            hashes[str(path)] = actual
        rst = Rst7.open(str(restart_path))
        xyz = np.asarray(rst.coordinates).reshape(-1, 3)
        force = np.loadtxt(force_path, skiprows=2, usecols=(1, 2, 3))
        if force.shape != (4662, 3):
            raise ValueError("Unexpected native force shape")
        bias = evaluate_bias(xyz, config)
        physical_gradient = force * SCALE
        gradient = physical_gradient + bias.gradient[active]
        energy = record["energy_hartree"] * EH + bias.energy
        imax = int(np.argmax(np.abs(gradient)) // 3)
        atom = model.topology.atoms[active[imax]]
        cutoff = cutoff_state(xyz, np.asarray(rst.box), model.qm)
        row = {
            "evaluation": number,
            "energy_hartree": record["energy_hartree"],
            "bias_kcal_mol": bias.energy,
            "total_kcal_mol": energy,
            "q": bias.q,
            "r_A": bias.r,
            "max_gradient_component_kcal_mol_A": float(np.abs(gradient).max()),
            "max_atom_1based": atom.idx + 1,
            "max_atom_name": atom.name,
            "max_atom_residue": f"{atom.residue.name}{atom.residue.idx + 1}",
            "qmcut_members": len(cutoff["members"]),
        }
        state = dict(
            row,
            xyz=xyz[active].copy(),
            gradient=gradient.copy(),
            physical_gradient=physical_gradient.copy(),
            cutoff=cutoff,
        )
        if previous is not None:
            delta = state["xyz"] - previous["xyz"]
            work = float(np.sum((gradient + previous["gradient"]) * delta) / 2)
            physical_work = float(
                np.sum((physical_gradient + previous["physical_gradient"]) * delta) / 2
            )
            de = energy - previous["total_kcal_mol"]
            row.update(
                {
                    "previous_evaluation": previous["evaluation"],
                    "max_component_displacement_A": float(np.abs(delta).max()),
                    "displacement_norm_A": float(np.linalg.norm(delta)),
                    "dE_kcal_mol": de,
                    "trapezoid_work_kcal_mol": work,
                    "dE_minus_work_kcal_mol": de - work,
                    "physical_dE_kcal_mol": (
                        record["energy_hartree"] - previous["energy_hartree"]
                    )
                    * EH,
                    "physical_trapezoid_work_kcal_mol": physical_work,
                    "gradient_change_max_component": float(
                        np.abs(gradient - previous["gradient"]).max()
                    ),
                }
            )
            changed = sorted(cutoff["members"] ^ previous["cutoff"]["members"])
            row["qmcut_membership_changes"] = len(changed)
            if changed:
                cutoff_events.append(
                    {
                        "from": previous["evaluation"],
                        "to": number,
                        "atoms": [
                            {
                                "id_1based": i + 1,
                                "name": model.topology.atoms[i].name,
                                "residue": model.topology.atoms[i].residue.name,
                                "residue_1based": model.topology.atoms[i].residue.idx
                                + 1,
                                "before_distance_A": previous["cutoff"]["boundary"].get(
                                    i
                                ),
                                "after_distance_A": cutoff["boundary"].get(i),
                                "entered": i in cutoff["members"],
                            }
                            for i in changed
                        ],
                    }
                )
        rows.append(row)
        states.append(state)
        previous = state
    comparisons = []
    last = states[-1]
    for state in states[-20:-1]:
        delta = last["xyz"] - state["xyz"]
        work = float(np.sum((last["gradient"] + state["gradient"]) * delta) / 2)
        comparisons.append(
            {
                "from": state["evaluation"],
                "to": last["evaluation"],
                "max_component_displacement_A": float(np.abs(delta).max()),
                "dE_kcal_mol": last["total_kcal_mol"] - state["total_kcal_mol"],
                "work_kcal_mol": work,
                "residual_kcal_mol": last["total_kcal_mol"]
                - state["total_kcal_mol"]
                - work,
                "gradient_change_max_component": float(
                    np.abs(last["gradient"] - state["gradient"]).max()
                ),
                "qmcut_changed_atoms_1based": [
                    i + 1
                    for i in sorted(
                        last["cutoff"]["members"] ^ state["cutoff"]["members"]
                    )
                ],
            }
        )
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (OUT / "steps.csv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "scope": "Retained data only; no new QM/MM calls or finite differences",
        "conversion": {"Eh_to_kcal_mol": EH, "raw_force_to_gradient": SCALE},
        "cutoff_definition": "MM atom within 10 A minimum-image distance of any real QM atom; ignores link-atom and residue policy",
        "limitations": [
            "Native forces printed to 7 decimal places; energy to 7 decimal places in Eh",
            "Cutoff geometric correlations alone do not establish engine implementation causality",
        ],
        "source_and_input_sha256": hashes,
        "steps": rows,
        "cutoff_events": cutoff_events,
        "comparisons_to_last": comparisons,
    }
    with (OUT / "audit.json").open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    logging.info(
        "Audited %d retained calls; %d cutoff events", len(rows), len(cutoff_events)
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
