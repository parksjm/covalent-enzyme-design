"""Compare the retained force jump with the identified MM cutoff pair."""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from itertools import product
from pathlib import Path

import numpy as np
from parmed import unit
from parmed.amber import Rst7
from parmed.geometry import box_lengths_and_angles_to_vectors

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = ROOT / (
    "qmmm/pa6_flipped_model30/07_ash274_comparison/"
    "08_activation_attack_grid/occupancy_callback_pilot"
)
sys.path[:0] = [str(ROOT / "src"), str(ROOT), str(BASE)]
from runtime_support import load_model  # noqa: E402


def main() -> None:
    """Write pair energy and Cartesian-gradient removal predictions."""
    model = load_model()
    out = BASE / "line_search_diagnosis"
    calls = BASE / "pilot/point_A"
    a, b = 11271, 11329
    rst = Rst7.open(str(calls / "eval_00283.input.ncrst"))
    xyz = np.asarray(rst.coordinates).reshape(-1, 3)
    cell = np.asarray(
        box_lengths_and_angles_to_vectors(*rst.box).value_in_unit(unit.angstroms)
    )
    shifts = np.asarray(list(product((-1, 0, 1), repeat=3))) @ cell
    deltas = xyz[a] - xyz[b] + shifts
    delta = deltas[np.argmin(np.linalg.norm(deltas, axis=1))]
    r = float(np.linalg.norm(delta))
    text = (calls / "eval_00284.mdout").read_text(errors="replace")
    alpha = float(re.search(r"Ewald Coefficient\s*=\s*([\d.]+)", text)[1])
    top = model.topology
    first, second = top.atoms[a], top.atoms[b]
    pair = (
        top.parm_data["NONBONDED_PARM_INDEX"][
            top.ptr("ntypes") * (first.nb_idx - 1) + second.nb_idx - 1
        ]
        - 1
    )
    acoef = top.parm_data["LENNARD_JONES_ACOEF"][pair]
    bcoef = top.parm_data["LENNARD_JONES_BCOEF"][pair]
    charge = 332.05221729 * first.charge * second.charge
    coulomb = charge * math.erfc(alpha * r) / r
    lj = acoef / r**12 - bcoef / r**6
    slope = charge * (
        -2 * alpha / math.sqrt(math.pi) * math.exp(-((alpha * r) ** 2)) / r
        - math.erfc(alpha * r) / r**2
    )
    slope += -12 * acoef / r**13 + 6 * bcoef / r**7
    predicted_a_change = -slope * delta / r
    forces = [
        np.loadtxt(
            calls / f"eval_{i:05d}.path_force.xyz", skiprows=2, usecols=(1, 2, 3)
        )
        for i in (283, 284)
    ]
    change = (forces[1] - forces[0]) * (-627.5094740631 / 0.529177**2)
    active = np.flatnonzero(model.active)
    maxrows = np.argsort(np.max(np.abs(change), axis=1))[-5:][::-1]
    rows = [
        {
            "atom_1based": int(active[i] + 1),
            "gradient_change_kcal_mol_A": change[i].tolist(),
        }
        for i in maxrows
    ]
    observed = [change[np.flatnonzero(active == atom)[0]].tolist() for atom in (a, b)]
    result = {
        "pair_atoms_1based": [a + 1, b + 1],
        "distance_A": r,
        "alpha_inverse_A_printed": alpha,
        "pair_coulomb_kcal_mol": coulomb,
        "pair_LJ_kcal_mol": lj,
        "predicted_removal_dE_kcal_mol": -coulomb - lj,
        "observed_dE_kcal_mol": -0.0010040141642093658,
        "native_energy_output_increment_kcal_mol": 1e-7 * 627.5094740631,
        "raw_force_component_output_increment_kcal_mol_A": 1e-7
        * 627.5094740631
        / 0.529177**2,
        "predicted_gradient_change_pair": [
            predicted_a_change.tolist(),
            (-predicted_a_change).tolist(),
        ],
        "observed_gradient_change_pair": observed,
        "top_gradient_changes": rows,
        "limitations": [
            "Alpha printed to5decimals; native forces to7decimals; native energies to7decimals Eh",
            "Standard erfc-screened Coulomb plus unshifted LJ evaluated from actual topology; engine not rerun",
        ],
        "final_analysis_source_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob("*.py")
        },
    }
    with (out / "pair_energy.json").open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
