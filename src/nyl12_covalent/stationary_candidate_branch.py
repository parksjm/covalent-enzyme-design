"""Diagnostic native descents without asserting a verified transition state."""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess

import numpy as np
from parmed.amber import Rst7

from nyl12_covalent.dftb_scc_validation import audit_dftb_scc
from nyl12_covalent.stationary_endpoints import classify_distances
from nyl12_covalent.stationary_native import NativeState, digest, write_json
from scripts.prepare_dimer_continuation import write_restart


def candidate_branch(args: tuple[Path, Path, str, float]) -> dict[str, object]:
    """Run one independent candidate-connectivity descent and retain raw evidence.

    Args:
        args: Frozen source campaign, new output, branch label and center energy.

    Returns:
        Separate process, force, energy and geometric diagnostic results.
    """
    source, output, label, center_energy = args
    os.nice(max(0, 10 - os.getpriority(os.PRIO_PROCESS, 0)))
    branch = output / "candidate_descents" / label
    branch.mkdir(parents=True)
    state = NativeState(source)
    sign = 1 if label == "plus" else -1
    mode = np.load(source / "ts_hessian/lowest_eigenvectors.npy")[:, 0]
    xyz = state.xyz.copy()
    xyz[state.active] += sign * 0.05 * (mode * state.inv_sqrt).reshape(-1, 3)
    write_restart(branch / "input.ncrst", xyz, state.restart)
    mdin = re.sub(r"\bmaxcyc\s*=\s*\d+", "maxcyc=20000", state.mdin)
    mdin = re.sub(r"\bmaxcycle\s*=\s*\d+", "maxcycle=20000", mdin)
    mdin = re.sub(r"\btol\s*=\s*[^,\s]+", "tol=1.d-5", mdin)
    mdin = re.sub(r"\btole\s*=\s*[^,\s]+", "tole=1.d-7", mdin)
    if (
        "optalg='LBFGS'" not in mdin
        or "dimer=.false." not in mdin
        or "dftb_telec=0.0" not in mdin
        or len(state.indices) != 4165
        or "qmcharge=-2" not in mdin.replace(" ", "")
    ):
        raise ValueError(
            "Candidate descent must retain declared state and minimum LBFGS"
        )
    (branch / "minimize.mdin").write_text(mdin)
    (branch / "topology.prmtop").symlink_to(state.topology)
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = (
        f"{source / 'library_build'}:/data/jpc/code/ambertools26/venv/lib"
    )
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
    write_json(
        branch / "specification.json",
        {
            "command": command,
            "mode": "candidate_connectivity_diagnostic",
            "verified_ts": False,
            "formal_IRC": False,
            "full_hessian_gate_passed": False,
            "source": str(source),
            "mode_sha256": digest(source / "ts_hessian/lowest_eigenvectors.npy"),
            "signed_initial_displacement": sign * 0.05,
            "maxcycle": 20000,
            "resume": "Preserve this attempt; use optimized.ncrst only after inspection, in a new directory.",
            "trajectory": "All native .path.xyz/.path_force.xyz retained without deletion.",
        },
    )
    try:
        with (branch / "native.log").open("w") as stream:
            process = subprocess.run(
                command,
                cwd=branch,
                env=env,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
        text = (branch / "minimize.mdout").read_text()
        electronic_history = audit_dftb_scc(text)
        electronic = bool(electronic_history["passed"])
        result: dict[str, object] = {
            "process_exit_status": process.returncode,
            "electronic_gate": electronic,
            "electronic_history": electronic_history,
            "native_optimizer_converged": "Optimisation converged" in text,
            "verified_ts": False,
            "formal_IRC": False,
        }
        if (
            process.returncode == 0
            and electronic
            and (branch / "optimized.ncrst").exists()
        ):
            final = np.asarray(
                Rst7.open(str(branch / "optimized.ncrst")).coordinates
            ).reshape(-1, 3)
            energy, gradient = state.gradient(final, branch / "final_gradient")
            distances = {
                name: float(np.linalg.norm(final[i] - final[j]))
                for name, (i, j) in {
                    "OG_C": (13915, 19915),
                    "OG_H": (13915, 13916),
                    "N_H": (19928, 13916),
                    "C_N": (19915, 19928),
                }.items()
            }
            maximum = float(np.max(np.abs(gradient)) * 0.529177)
            result.update(
                {
                    "energy_hartree": energy,
                    "energy_below_candidate": energy < center_energy,
                    "max_gradient_hartree_per_bohr": maximum,
                    "force_converged": maximum <= 1e-5,
                    "distances_angstrom": distances,
                    "geometry_class": classify_distances(distances),
                }
            )
        write_json(branch / "diagnostic_acceptance.json", result)
        return result
    except Exception as exc:
        write_json(branch / "failure.json", {"error": str(exc), "verified_ts": False})
        raise
