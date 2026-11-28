"""Strict electronic and divergence guards for bounded relay native searches."""

from __future__ import annotations

from datetime import datetime, timezone
import os
from pathlib import Path
import re
import subprocess
import time
from typing import Any

import numpy as np
from parmed.amber import Rst7

from .dftb_scc_validation import audit_dftb_scc
from .nterminal_relay_native import hydrogen_audit, unchanged_exterior
from .stationary_native import NativeState, digest, write_json
from .relay_release import geometry
from .relay_scc_guard import scc_guard_reason


def invalid_trial(
    text: str, initial_energy: float, allowance: float, *, final: bool = False
) -> str | None:
    """Return the first strict electronic or physical-energy divergence reason.

    Args:
        text: Current native output, possibly from an active process.
        initial_energy: Independently evaluated initial energy in Hartree.
        allowance: Maximum native total energy rise in kcal/mol.
        final: Reject unresolved warnings when no further recovery can occur.

    Returns:
        Explicit failure reason, or None when no guard has fired.
    """
    reason = scc_guard_reason(text, final=final, maximum_iterations=500)
    if reason is not None:
        return reason
    energies = re.findall(r"Etot\s*=\s*([-+\d.EeDd]+)", text)
    if energies:
        values = np.array(
            [float(v.replace("D", "E").replace("d", "e")) for v in energies]
        )
        if (
            not np.isfinite(values).all()
            or np.max(values) > initial_energy * 627.5094740631 + allowance
        ):
            return "Native physical energy exceeded declared initial-energy allowance"
    return None


def run_guarded(
    model: Path,
    start: Path,
    output: Path,
    mdin: str,
    *,
    allowance: float = 25.0,
    hours: float = 4.0,
) -> dict[str, Any]:
    """Run an immutable native stage with bounded SCC recovery and valid forces.

    Args:
        model: Original unchanged QM/MM model and precision library.
        start: Independently SCC-screened original full-system seed.
        output: Fresh native attempt directory.
        mdin: Explicit native dimer or minimization input.
        allowance: Allowed energy rise in kcal/mol before stopping exact child.
        hours: Hard process wall-time budget.

    Returns:
        Separate process, force, electronic, exterior and chemical acceptance.
    """
    if allowance <= 0 or hours <= 0:
        raise ValueError("Positive energy and process budgets required")
    output.mkdir(parents=True, exist_ok=False)
    state = NativeState(model, start)
    reference = NativeState(model)
    unchanged_exterior(reference.xyz, state.xyz, state.active)
    for key in ("active", "qmmask", "qm_theory"):
        pattern = rf"\b{key}\s*=\s*'([^']+)'"
        if re.findall(pattern, mdin, re.I) != re.findall(
            pattern, state.mdin_source, re.I
        ):
            raise ValueError(f"Changed original {key}")
    for key in (
        "dftb_telec",
        "qmcharge",
        "scfconv",
        "cut",
        "qmcut",
        "qm_ewald",
        "qm_pme",
        "qmmm_int",
        "adjust_q",
        "tight_p_conv",
    ):
        pattern = rf"\b{key}\s*=\s*([-+\d.EeDd]+)"
        if re.findall(pattern, mdin, re.I) != re.findall(
            pattern, state.mdin_source, re.I
        ):
            raise ValueError(f"Changed original Hamiltonian control {key}")
    for key in ("ntr", "nmropt"):
        if re.findall(rf"\b{key}\s*=\s*(\d+)", mdin, re.I) != ["0"]:
            raise ValueError(f"Unexpected {key} restraint")
    initial_energy, initial_gradient = state.gradient(
        state.xyz, output / "initial_gradient"
    )
    (output / "input.mdin").write_text(mdin)
    command = [
        str(state.sander),
        "-O",
        "-i",
        "input.mdin",
        "-p",
        str(state.topology),
        "-c",
        str(start),
        "-o",
        "native.mdout",
        "-r",
        "optimized.ncrst",
        "-x",
        "optimization.nc",
        "-inf",
        "mdinfo",
    ]
    env = dict(
        os.environ,
        AMBERHOME="/data/jpc/code/ambertools26/venv",
        LD_LIBRARY_PATH=f"{model / 'library_build'}:/data/jpc/code/ambertools26/venv/lib",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    process = {
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "initial_energy_hartree": initial_energy,
        "initial_max_gradient_Eh_per_bohr": float(
            np.max(np.abs(initial_gradient)) * 0.529177
        ),
        "initial_geometry": geometry(state.xyz),
        "energy_rise_allowance_kcal_mol": allowance,
        "wall_hours": hours,
        "scc_recovery_policy": {
            "maximum_iterations_per_evaluation": 500,
            "allow_native_same_evaluation_recovery": True,
            "returned_and_final_gradients_must_pass_exact_scc_audit": True,
            "physical_SCF_tolerances_and_temperature_unchanged": True,
        },
        "input_sha256": {
            str(p): digest(p)
            for p in (
                start,
                state.topology,
                output / "input.mdin",
                state.sander,
                model / "library_build/libdlfind.so",
            )
        },
    }
    write_json(output / "process.json", process)
    started = time.monotonic()
    reason = None
    with (output / "engine.log").open("x") as stream:
        child = subprocess.Popen(
            command, cwd=output, env=env, stdout=stream, stderr=subprocess.STDOUT
        )
        process["pid"] = child.pid
        write_json(output / "process.json", process)
        while child.poll() is None:
            native = output / "native.mdout"
            if native.exists():
                reason = invalid_trial(
                    native.read_text(errors="replace"), initial_energy, allowance
                )
            if time.monotonic() - started > hours * 3600:
                reason = "Declared hard wall-time budget exceeded"
            if reason is not None:
                child.terminate()
                try:
                    child.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
                break
            time.sleep(5)
    process.update(
        exit_status=child.returncode,
        abort_reason=reason,
        ended_utc=datetime.now(timezone.utc).isoformat(),
    )
    write_json(output / "process.json", process)
    text = (output / "native.mdout").read_text(errors="replace")
    reason = reason or invalid_trial(text, initial_energy, allowance, final=True)
    result = {
        "process_completed": child.returncode == 0,
        "abort_reason": reason,
        "native_converged": "Optimisation converged" in text,
        "full_active_index_verified": False,
        "accepted_force_candidate": False,
    }
    write_json(output / "scientific_acceptance.json", result)
    if child.returncode or reason or not result["native_converged"]:
        raise RuntimeError(f"Guarded native stage failed: {result}")
    xyz = np.asarray(Rst7.open(str(output / "optimized.ncrst")).coordinates).reshape(
        (-1, 3)
    )
    exterior = unchanged_exterior(reference.xyz, xyz, state.active)
    energy, gradient = state.gradient(xyz, output / "final_gradient")
    maximum = float(np.max(np.abs(gradient)) * 0.529177)
    hydrogens = hydrogen_audit(model, xyz)
    dimer = bool(re.search(r"dimer\s*=\s*\.true\.", mdin, re.I))
    accepted = maximum <= (1e-4 if dimer else 1e-5) and hydrogens["passed"]
    result.update(
        accepted_force_candidate=accepted,
        energy_hartree=energy,
        maximum_gradient_Eh_per_bohr=maximum,
        frozen_max_change_A=exterior,
        geometry=geometry(xyz),
        nonrelay_hydrogens=hydrogens,
        electronic_history=audit_dftb_scc(text),
    )
    write_json(output / "scientific_acceptance.json", result)
    if not accepted:
        raise RuntimeError("Independent final gradient or nonrelay proton gate failed")
    return result
