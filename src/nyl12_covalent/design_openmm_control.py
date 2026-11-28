"""Run a bounded, unrestrained OpenMM minimization of a prepared E-S model."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
from pathlib import Path
from typing import Any

FORCE_TOLERANCE_KJ_MOL_NM = 2.092
MAX_ITERATIONS = 10_000
PLATFORM_PRECISION = "mixed"
PME_CUTOFF_NM = 1.0
LOGGER = logging.getLogger(__name__)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_energy_and_rms_force(state: Any, unit: Any) -> tuple[float, float, int]:
    """Return finite energy and all-components RMS force in OpenMM units."""
    import numpy as np

    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    if forces.ndim != 2 or forces.shape[1] != 3 or not forces.size:
        raise ValueError("OpenMM returned forces with an invalid shape")
    if not math.isfinite(energy) or not np.isfinite(forces).all():
        raise ValueError("OpenMM returned non-finite energy or forces")
    rms_force = float(np.sqrt(np.mean(forces * forces)))
    if not math.isfinite(rms_force):
        raise ValueError("OpenMM returned a non-finite RMS force")
    return energy, rms_force, int(forces.shape[0])


def run_minimization(topology: Path, coordinates: Path, output: Path) -> dict[str, Any]:
    """Minimize an Amber topology/restart on visible CUDA device zero.

    Args:
        topology: Prepared Amber ``prmtop`` with validated E-S parameters.
        coordinates: Matching Amber ``rst7`` containing a periodic box.
        output: New directory for coordinates and the acceptance receipt.

    Returns:
        Receipt separating process completion, energy change, and force
        convergence. Coordinates are saved in nanometers in NumPy format.

    Raises:
        FileExistsError: If the output directory already exists.
        FileNotFoundError: If an input file is missing.
        ValueError: If the GPU mapping, periodic box, or numerical results fail.
    """
    if not topology.is_file():
        raise FileNotFoundError(topology)
    if not coordinates.is_file():
        raise FileNotFoundError(coordinates)
    if not topology.is_absolute() or not coordinates.is_absolute():
        raise ValueError("Topology and coordinates must be absolute paths")
    if not output.is_absolute():
        raise ValueError("Output must be an absolute path")
    if output.exists():
        raise FileExistsError(output)
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "0":
        raise ValueError("CUDA_VISIBLE_DEVICES must be exactly '0' for this control")

    import numpy as np
    import openmm
    from openmm import app, unit

    LOGGER.info("Starting OpenMM minimization for %s", topology)
    amber_topology = app.AmberPrmtopFile(str(topology))
    amber_coordinates = app.AmberInpcrdFile(str(coordinates))
    if amber_coordinates.boxVectors is None:
        raise ValueError("Amber restart does not contain periodic box vectors")
    system = amber_topology.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=PME_CUTOFF_NM * unit.nanometer,
        constraints=None,
        rigidWater=False,
        ewaldErrorTolerance=5e-4,
    )
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName("CUDA")
    properties = {"DeviceIndex": "0", "Precision": PLATFORM_PRECISION}
    simulation = app.Simulation(
        amber_topology.topology, system, integrator, platform, properties
    )
    simulation.context.setPeriodicBoxVectors(*amber_coordinates.boxVectors)
    simulation.context.setPositions(amber_coordinates.positions)

    initial_state = simulation.context.getState(getEnergy=True, getForces=True)
    initial_energy, _, atom_count = _finite_energy_and_rms_force(initial_state, unit)
    openmm.LocalEnergyMinimizer.minimize(
        simulation.context,
        FORCE_TOLERANCE_KJ_MOL_NM * unit.kilojoule_per_mole / unit.nanometer,
        MAX_ITERATIONS,
    )
    final_state = simulation.context.getState(
        getEnergy=True, getForces=True, getPositions=True
    )
    final_energy, final_rms_force, final_atom_count = _finite_energy_and_rms_force(
        final_state, unit
    )
    positions_nm = np.asarray(
        final_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
        dtype=float,
    )
    if (
        positions_nm.ndim != 2
        or positions_nm.shape[1] != 3
        or positions_nm.shape[0] != final_atom_count
        or final_atom_count != atom_count
    ):
        raise ValueError("OpenMM returned positions with an invalid shape")
    if not np.isfinite(positions_nm).all():
        raise ValueError("OpenMM returned non-finite coordinates")

    output.mkdir(parents=True, exist_ok=False)
    np.save(output / "positions_nm.npy", positions_nm)
    with (output / "final.pdb").open("w", encoding="utf-8") as stream:
        app.PDBFile.writeFile(
            amber_topology.topology, final_state.getPositions(), stream, keepIds=True
        )
    energy_change = final_energy - initial_energy
    receipt: dict[str, Any] = {
        "status": "minimization_complete",
        "topology": str(topology),
        "topology_sha256": _sha256(topology),
        "coordinates": str(coordinates),
        "coordinates_sha256": _sha256(coordinates),
        "output_positions": str(output / "positions_nm.npy"),
        "output_pdb": str(output / "final.pdb"),
        "platform": "CUDA",
        "device_index": 0,
        "cuda_visible_devices": "0",
        "precision": PLATFORM_PRECISION,
        "nonbonded_method": "PME",
        "cutoff_nm": PME_CUTOFF_NM,
        "constraints": None,
        "rigid_water": False,
        "restraints": "none",
        "dynamics_performed": False,
        "max_iterations": MAX_ITERATIONS,
        "force_tolerance_kj_mol_nm": FORCE_TOLERANCE_KJ_MOL_NM,
        "initial_energy_kj_mol": initial_energy,
        "final_energy_kj_mol": final_energy,
        "energy_change_kj_mol": energy_change,
        "energy_decreased": energy_change < 0.0,
        "final_rms_force_kj_mol_nm": final_rms_force,
        "final_rms_gradient_kcal_mol_angstrom": final_rms_force / 41.84,
        "force_converged": final_rms_force <= FORCE_TOLERANCE_KJ_MOL_NM,
        "coordinates_finite": True,
    }
    receipt_path = output / "acceptance.json"
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    LOGGER.info(
        "Finished OpenMM minimization; force-converged=%s, RMS=%.6g kcal/mol/A",
        receipt["force_converged"],
        receipt["final_rms_gradient_kcal_mol_angstrom"],
    )
    return receipt


def main() -> None:
    """Run one explicitly named OpenMM control minimization."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run_minimization(args.topology, args.coordinates, args.output)


if __name__ == "__main__":
    main()
