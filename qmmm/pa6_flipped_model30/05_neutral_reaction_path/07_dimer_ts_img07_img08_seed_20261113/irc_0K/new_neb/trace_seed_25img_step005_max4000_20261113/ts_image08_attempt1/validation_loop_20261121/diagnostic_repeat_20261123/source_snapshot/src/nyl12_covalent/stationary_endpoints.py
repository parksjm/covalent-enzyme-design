"""Both-sign downhill endpoint minimization from a verified Hessian mode."""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess

import numpy as np
from parmed.amber import Rst7

from nyl12_covalent.dftb_scc_validation import audit_dftb_scc
from nyl12_covalent.stationary_native import NativeState, write_json
from scripts.prepare_dimer_continuation import write_restart


def classify_distances(distances: dict[str, float]) -> str:
    """Classify the declared direct-route bond/proton inventory geometrically."""
    if distances["C_N"] < 1.6 and distances["OG_C"] > 2.3 and distances["OG_H"] < 1.2:
        return "reactant_ES"
    if distances["C_N"] > 2.0 and distances["OG_C"] < 1.7 and distances["N_H"] < 1.2:
        return "acyl_enzyme"
    return "unassigned"


def endpoints(root: Path, ts_energy: float) -> list[Path]:
    """Minimize both directions along the actual verified negative eigenvector.

    Args:
        root: Validation campaign with an accepted complete TS Hessian.
        ts_energy: Native stationary-candidate energy in Hartree.

    Returns:
        Two force-converged endpoint restarts if chemical connectivity passes.
    """
    state = NativeState(root)
    modes = np.load(root / "ts_hessian/lowest_eigenvectors.npy")
    direction = (modes[:, 0] * state.inv_sqrt).reshape((-1, 3))
    out = root / "endpoints"
    out.mkdir()
    points = []
    labels = []
    for sign, label in ((-1, "minus"), (1, "plus")):
        branch = out / label
        branch.mkdir()
        accepted = None
        for displacement in (0.05, 0.1, 0.2, 0.4):
            xyz = state.xyz.copy()
            xyz[state.active] += sign * displacement * direction
            energy, _ = state.gradient(xyz, branch / f"displace_{displacement}")
            if energy < ts_energy - 1e-7:
                accepted = xyz
                break
        if accepted is None:
            raise RuntimeError(f"{label} branch has no demonstrated downhill start")
        write_restart(branch / "input.ncrst", accepted, state.restart)
        mdin = re.sub(r"maxcyc=1\b", "maxcyc=20000", state.mdin)
        mdin = mdin.replace("maxcycle=1,", "maxcycle=20000,")
        mdin = mdin.replace("tol=1.d-4, tole=1.d-5", "tol=1.d-5, tole=1.d-7")
        (branch / "minimize.mdin").write_text(mdin)
        env = dict(os.environ)
        env["LD_LIBRARY_PATH"] = (
            f"{root / 'library_build'}:/data/jpc/code/ambertools26/venv/lib"
        )
        (branch / "topology.prmtop").symlink_to(state.topology)
        command = [
            str(state.sander),
            "-O",
            "-i",
            "minimize.mdin",
            "-p",
            "topology.prmtop",
            "-c",
            "input.ncrst",
            "-o",
            "minimize.mdout",
            "-r",
            "optimized.ncrst",
            "-inf",
            "mdinfo",
        ]
        with (branch / "native.log").open("w") as log:
            subprocess.run(
                command,
                cwd=branch,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        text = (branch / "minimize.mdout").read_text()
        electronic_history = audit_dftb_scc(text)
        if "Optimisation converged" not in text or not electronic_history["passed"]:
            raise RuntimeError(f"{label}: native force/electronic convergence failed")
        restart = branch / "optimized.ncrst"
        xyz = np.asarray(Rst7.open(str(restart)).coordinates).reshape((-1, 3))
        energy, gradient = state.gradient(xyz, branch / "final_gradient")
        max_gradient = float(np.max(np.abs(gradient)) * 0.529177)
        if max_gradient > 1e-5 or energy >= ts_energy:
            raise RuntimeError(f"{label}: endpoint force/energy acceptance failed")
        pairs = {
            "OG_C": (13915, 19915),
            "OG_H": (13915, 13916),
            "N_H": (19928, 13916),
            "C_N": (19915, 19928),
        }
        distances = {
            name: float(np.linalg.norm(xyz[i] - xyz[j]))
            for name, (i, j) in pairs.items()
        }
        chemistry = classify_distances(distances)
        write_json(
            branch / "scientific_acceptance.json",
            {
                "force_converged": True,
                "energy_hartree": energy,
                "max_gradient_hartree_per_bohr": max_gradient,
                "distances_angstrom": distances,
                "chemistry": chemistry,
                "negative_mode_initial_displacement_mass_weighted": sign * displacement,
                "description": "Two-sided downhill eigenmode displacements followed by native LBFGS; local connectivity, not formal mass-weighted IRC.",
                "electronic_history": electronic_history,
            },
        )
        points.append(restart)
        labels.append(chemistry)
    if set(labels) != {"reactant_ES", "acyl_enzyme"}:
        raise RuntimeError(f"Endpoints do not establish intended chemistry: {labels}")
    write_json(
        out / "connectivity_acceptance.json",
        {
            "passed": True,
            "endpoint_labels": labels,
            "source_ts": str(root / "inputs/optimized.ncrst"),
            "mode_source": str(root / "ts_hessian/lowest_eigenvectors.npy"),
        },
    )
    return points
