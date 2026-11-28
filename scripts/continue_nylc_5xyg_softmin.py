#!/usr/bin/env python3
"""Continue the failed 5XYG repair with frozen original atoms and native soft potential."""

from __future__ import annotations

import argparse
import json
import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
RUN = (
    ROOT
    / "analysis/negoro_2023_protonation/pdbfixer_repair_20261120/repair_run_20261121"
)
LOGGER = logging.getLogger(__name__)


def expected_sequences(path: Path) -> dict[str, list[str]]:
    """Read split-chain SEQRES records from the retained PDBFixer input."""
    sequences: dict[str, list[str]] = defaultdict(list)
    for line in path.read_text().splitlines():
        if line.startswith("SEQRES"):
            sequences[line[11:12]].extend(line[19:70].split())
    if set(sequences) != set("ABCD"):
        raise ValueError(f"Expected four split SEQRES chains, got {sorted(sequences)}")
    return dict(sequences)


def continue_softmin(
    candidate: Path, original_input: Path, destination: Path
) -> dict[str, Any]:
    """Minimize only newly modeled heavy atoms using PDBFixer's native soft field."""
    import openmm as mm
    from openmm import Platform, unit
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    from repair_nylc_5xyg_pdbfixer import read_atoms, validate_repair

    destination.parent.mkdir(parents=True, exist_ok=True)
    original = read_atoms(original_input)
    fixer = PDBFixer(filename=str(candidate))
    fixer.platform = Platform.getPlatformByName("CPU")
    topology = fixer.topology
    atom_keys = {
        (
            atom.residue.chain.id,
            int(atom.residue.id),
            atom.residue.insertionCode or " ",
            atom.name,
        ): atom
        for atom in topology.atoms()
    }
    missing = set(original) - set(atom_keys)
    if missing:
        raise ValueError(f"Candidate omitted original atoms: {sorted(missing)[:5]}")
    new_atoms = [atom for key, atom in atom_keys.items() if key not in original]
    if not new_atoms:
        raise ValueError("Candidate contains no newly modeled atoms to optimize")

    forcefield = fixer._createForceField(topology, False)
    system = forcefield.createSystem(topology)
    for key in original:
        system.setParticleMass(atom_keys[key].index, 0 * unit.dalton)
    nonbonded = [
        force
        for force in system.getForces()
        if isinstance(force, mm.CustomNonbondedForce)
    ]
    if len(nonbonded) != 1:
        raise ValueError(
            f"Expected one native soft CustomNonbondedForce, got {len(nonbonded)}"
        )
    nonbonded[0].addInteractionGroup(
        {atom.index for atom in new_atoms}, set(range(system.getNumParticles()))
    )
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(system, integrator, fixer.platform)
    context.setPositions(fixer.positions)
    initial_energy = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoules_per_mole)
    )
    if not math.isfinite(initial_energy):
        raise ValueError(
            f"Initial soft-potential energy is nonfinite: {initial_energy}"
        )
    mm.LocalEnergyMinimizer.minimize(
        context,
        tolerance=0.1 * unit.kilojoules_per_mole / unit.nanometer,
        maxIterations=20000,
    )
    state = context.getState(getPositions=True, getEnergy=True, getForces=True)
    final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    forces = state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoules_per_mole / unit.nanometer
    )
    if not math.isfinite(final_energy) or not all(
        math.isfinite(float(value)) for value in forces.flat
    ):
        raise ValueError("Soft minimization produced nonfinite energy or forces")
    fixer.positions = state.getPositions()
    with destination.open("w") as handle:
        PDBFile.writeFile(topology, fixer.positions, handle, keepIds=True)

    repaired = PDBFile(str(destination))
    validation = validate_repair(
        original,
        destination,
        expected_sequences(original_input),
        repaired.topology,
    )
    report = {
        "state": "accepted"
        if validation["maximum_original_heavy_displacement_angstrom"] <= 0.002
        else "rejected",
        "input_candidate": str(candidate),
        "original_split_input": str(original_input),
        "output": str(destination),
        "platform": fixer.platform.getName(),
        "forcefield": "PDBFixer v1.12 native soft.xml via PDBFixer._createForceField(topology, False)",
        "new_atom_count": len(new_atoms),
        "frozen_original_atom_count": len(original),
        "frozen_original_mass_dalton": 0,
        "nonbonded_force_count": len(nonbonded),
        "minimizer": {
            "name": "OpenMM LocalEnergyMinimizer",
            "tolerance_kj_mol_nm": 0.1,
            "max_iterations": 20000,
        },
        "initial_soft_energy_kj_mol": initial_energy,
        "final_soft_energy_kj_mol": final_energy,
        "max_force_kj_mol_nm": float(
            max(math.sqrt(sum(float(x) ** 2 for x in row)) for row in forces)
        ),
        "energy_interpretation": "PDBFixer soft-potential diagnostic values, not physical protein MM energies",
        **validation,
    }
    report_path = destination.with_name("softmin_scientific_acceptance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info("Soft-potential continuation result: %s", report["state"])
    return report


def main() -> None:
    """Run the one permitted soft-potential continuation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--original-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    result = continue_softmin(args.candidate, args.original_input, args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
