"""Execute one native relay release and record independent scientific gates."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from parmed.amber import Rst7

from nyl12_covalent.dftb_scc_validation import audit_dftb_scc
from nyl12_covalent.relay_release import geometry
from nyl12_covalent.stationary_native import NativeState, digest, write_json


def run(seed: Path) -> None:
    """Run a fresh unbiased native minimization and retain process/chemistry gates.

    Args:
        seed: Prepared immutable per-seed root within the bounded campaign.
    """
    state = NativeState(seed)
    manifest = json.loads((seed.parent / "prepared.json").read_text())
    for name, expected in manifest["parameter_sha256"].items():
        installed = Path(os.environ["AMBERHOME"]) / "dat/slko/3ob-3-1" / name
        if digest(installed) != expected:
            raise ValueError("Installed DFTB parameter differs from prepared snapshot")
    initial_energy, initial_gradient = state.gradient(
        state.xyz, seed / "initial_gradient", precision=False
    )
    write_json(
        seed / "initial_acceptance.json",
        {
            "energy_hartree": initial_energy,
            "max_gradient_hartree_per_bohr": float(
                np.max(np.abs(initial_gradient)) * 0.529177
            ),
            "geometry": geometry(state.xyz),
        },
    )
    command = [
        str(state.sander),
        "-O",
        "-i",
        "minimize.mdin",
        "-p",
        "inputs/system.prmtop",
        "-c",
        "inputs/optimized.ncrst",
        "-o",
        "minimize.mdout",
        "-r",
        "optimized.ncrst",
        "-x",
        "optimization.nc",
        "-inf",
        "mdinfo",
    ]
    start = time.time()
    with (seed / "native.log").open("x") as log:
        process = subprocess.run(
            command, cwd=seed, stdout=log, stderr=subprocess.STDOUT, check=False
        )
    write_json(
        seed / "process_acceptance.json",
        {
            "exit_status": process.returncode,
            "elapsed_seconds": time.time() - start,
            "command": command,
        },
    )
    text = (seed / "minimize.mdout").read_text()
    scc = audit_dftb_scc(text)
    native_converged = "Optimisation converged" in text
    if process.returncode or not scc["passed"] or not native_converged:
        write_json(
            seed / "scientific_acceptance.json",
            {
                "accepted": False,
                "native_converged": native_converged,
                "electronic_history": scc,
                "reason": "Process/SCC/native convergence gate failed; inspect retained trajectory",
            },
        )
        raise RuntimeError("Native optimization incomplete or electronic gate failed")
    final = np.asarray(Rst7.open(str(seed / "optimized.ncrst")).coordinates).reshape(
        (-1, 3)
    )
    energy, gradient = state.gradient(final, seed / "final_gradient", precision=False)
    maximum = float(np.max(np.abs(gradient)) * 0.529177)
    displacement = np.linalg.norm(final - state.xyz, axis=1)
    frozen_max = float(displacement[~state.active].max())
    result = {
        "accepted": maximum <= 1e-5 and energy < initial_energy and frozen_max < 1e-10,
        "force_converged": maximum <= 1e-5,
        "max_gradient_hartree_per_bohr": maximum,
        "native_converged": native_converged,
        "electronic_history": scc,
        "initial_energy_hartree": initial_energy,
        "final_energy_hartree": energy,
        "energy_decreased": energy < initial_energy,
        "frozen_max_displacement_A": frozen_max,
        "active_max_displacement_A": float(displacement[state.active].max()),
        "active_rms_displacement_A": float(
            np.sqrt(np.mean(displacement[state.active] ** 2))
        ),
        "initial_geometry": geometry(state.xyz),
        "final_geometry": geometry(final),
        "stationary_point_order": "Not evaluated: no Hessian performed",
        "chemical_retention": "Inspect final proton/bond geometry independently of numerical acceptance",
    }
    write_json(seed / "scientific_acceptance.json", result)
    if not result["accepted"]:
        raise RuntimeError("Independent force/energy/frozen-boundary gate failed")
