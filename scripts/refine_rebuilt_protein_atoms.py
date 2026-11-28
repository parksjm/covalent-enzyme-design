#!/usr/bin/env python3
"""Refine rebuilt protein atoms with Amber19SB/GBn2 and a frozen observed core."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atom_key(atom: Any) -> tuple[str, int, str, str]:
    """Return a persistent PDB atom identity."""
    return (
        atom.residue.chain.id,
        int(atom.residue.id),
        atom.residue.insertionCode or " ",
        atom.name,
    )


def assert_core_retained(original: dict, observed: dict, tolerance: float) -> float:
    """Require every experimental heavy atom within the stated displacement."""
    if set(original) - set(observed):
        raise ValueError("Experimental heavy atoms missing from model")
    maximum = max(math.dist(xyz, observed[key]) for key, xyz in original.items())
    if maximum > tolerance:
        raise ValueError(f"Experimental heavy atom moved {maximum:.6f} A")
    return maximum


def geometry_audit(
    path: Path, original: dict, expected: dict, termini: dict
) -> dict[str, Any]:
    """Check sequence, termini, peptide bonds, observed core and modeled clashes."""
    from openmm.app import PDBFile
    from postflight_nylc_5xyg_repair import clash_audit
    from repair_nylc_5xyg_pdbfixer import read_atoms

    pdb = PDBFile(str(path))
    observed = read_atoms(path)
    drift = assert_core_retained(original, observed, 0.002)
    actual = {c.id: [r.name for r in c.residues()] for c in pdb.topology.chains()}
    failures = []
    if actual != expected:
        failures.append("sequence_mismatch")
    peptides = []
    bond_keys = {frozenset((atom_key(a), atom_key(b))) for a, b in pdb.topology.bonds()}
    for chain in pdb.topology.chains():
        residues = list(chain.residues())
        declared = termini[chain.id]
        if [int(r.id) for r in residues] != list(
            range(declared["n"], declared["c"] + 1)
        ):
            failures.append(f"residue_inventory:{chain.id}")
        if (int(residues[0].id), int(residues[-1].id)) != (
            declared["n"],
            declared["c"],
        ):
            failures.append(f"terminus_mismatch:{chain.id}")
        if "N" not in {a.name for a in residues[0].atoms()}:
            failures.append(f"missing_N_terminal:{chain.id}")
        if sum(a.name == "OXT" for r in residues for a in r.atoms()) != 1:
            failures.append(f"OXT_inventory:{chain.id}")
        if "OXT" not in {a.name for a in residues[-1].atoms()}:
            failures.append(f"missing_OXT:{chain.id}")
        for left, right in zip(residues, residues[1:]):
            c = (chain.id, int(left.id), left.insertionCode or " ", "C")
            n = (chain.id, int(right.id), right.insertionCode or " ", "N")
            distance = math.dist(observed[c], observed[n])
            if not 1.1 <= distance <= 1.8 or frozenset((c, n)) not in bond_keys:
                peptides.append({"left": c, "right": n, "distance_A": distance})
    crosschain = [
        (atom_key(a), atom_key(b))
        for a, b in pdb.topology.bonds()
        if a.residue.chain != b.residue.chain and {a.name, b.name} == {"C", "N"}
    ]
    clashes = clash_audit(path, set(observed) - set(original))
    if peptides:
        failures.append("peptide_geometry_or_bond")
    if crosschain:
        failures.append("cross_chain_peptide")
    if clashes["clash_count"]:
        failures.append("modeled_atom_clashes")
    return {
        "failures": failures,
        "maximum_original_drift_A": drift,
        "bad_peptides": peptides,
        "crosschain_peptide_bonds": crosschain,
        "clash_audit": clashes,
        "sequences": actual,
    }


def refine(
    candidate: Path, original_path: Path, output_dir: Path, termini_path: Path
) -> dict[str, Any]:
    """Minimize rebuilt heavy atoms and temporary H while freezing observed atoms."""
    import openmm as mm
    from openmm import unit
    from openmm.app import ForceField, Modeller, NoCutoff, PDBFile
    from openmm.app import forcefield as ff_module
    from repair_nylc_5xyg_pdbfixer import read_atoms

    output_dir.mkdir(parents=True, exist_ok=False)
    report: dict[str, Any] = {
        "state": "running",
        "method": "Amber19SB/GBn2",
        "hydrogen_build_pH": 8.0,
        "hydrogen_interpretation": "temporary building state, not pKa microstate assignment",
        "openmm_version": mm.__version__,
        "cpu_threads": 2,
    }
    report_path = output_dir / "scientific_acceptance.json"
    try:
        for source in (
            candidate,
            original_path,
            Path(__file__),
            Path(__file__).with_name("repair_nylc_5xyg_pdbfixer.py"),
            Path(__file__).with_name("postflight_nylc_5xyg_repair.py"),
        ):
            shutil.copy2(source, output_dir / source.name)
        report["source_hashes"] = {
            str(p): sha256(p)
            for p in (
                candidate,
                original_path,
                Path(__file__),
                Path(__file__).with_name("repair_nylc_5xyg_pdbfixer.py"),
                Path(__file__).with_name("postflight_nylc_5xyg_repair.py"),
            )
        }
        data = Path(ff_module.__file__).parent / "data"
        ff_names = ["amber19-all.xml", "implicit/gbn2.xml"]
        forcefield = ForceField(*ff_names)
        # Preserve all installed XML parameter evidence, including transitive includes.
        report["forcefield_xml_hashes"] = {
            str(p.relative_to(data)): sha256(p) for p in sorted(data.rglob("*.xml"))
        }
        termini = json.loads(termini_path.read_text())
        report["declared_termini"] = termini
        report["termini_sha256"] = sha256(termini_path)
        shutil.copy2(termini_path, output_dir / "declared_termini.json")
        original = read_atoms(original_path)
        heavy = read_atoms(candidate)
        assert_core_retained(original, heavy, 0.002)
        pdb = PDBFile(str(candidate))
        expected: dict[str, list[str]] = defaultdict(list)
        for line in original_path.read_text().splitlines():
            if line.startswith("SEQRES"):
                expected[line[11:12]].extend(line[19:70].split())
        if not expected:
            raise ValueError(
                "Original input requires explicit SEQRES for sequence acceptance"
            )
        if set(termini) != set(expected):
            raise ValueError("Declared termini chains differ from SEQRES")
        modeller = Modeller(pdb.topology, pdb.positions)
        platform = mm.Platform.getPlatformByName("CPU")
        variants = modeller.addHydrogens(forcefield, pH=8.0, platform=platform)
        report["hydrogen_atom_inventory"] = [
            {
                "chain": r.chain.id,
                "residue": r.id,
                "hydrogens": [a.name for a in r.atoms() if a.element.symbol == "H"],
            }
            for r in modeller.topology.residues()
        ]
        report["actual_hydrogen_variants"] = [
            {"chain": r.chain.id, "residue": r.id, "name": r.name, "variant": v}
            for r, v in zip(modeller.topology.residues(), variants)
        ]
        positions = modeller.positions.value_in_unit(unit.angstrom)
        after_h = {
            atom_key(a): tuple(positions[a.index])
            for a in modeller.topology.atoms()
            if a.element.symbol != "H"
        }
        report["hydrogen_addition_maximum_heavy_drift_A"] = assert_core_retained(
            heavy, after_h, 1e-8
        )
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=NoCutoff,
            constraints=None,
            removeCMMotion=False,
        )
        movable = []
        mapping = []
        for atom in modeller.topology.atoms():
            key = atom_key(atom)
            frozen = key in original and atom.element.symbol != "H"
            if frozen:
                system.setParticleMass(atom.index, 0 * unit.dalton)
            else:
                movable.append(atom.index)
            mapping.append(
                {
                    "index": atom.index,
                    "identity": key,
                    "element": atom.element.symbol,
                    "frozen": frozen,
                }
            )
        (output_dir / "atom_mapping.json").write_text(
            json.dumps(mapping, indent=2) + "\n"
        )
        report["frozen_original_heavy_count"] = len(original)
        report["movable_count"] = len(movable)
        integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
        context = mm.Context(system, integrator, platform, {"Threads": "2"})
        context.setPositions(modeller.positions)
        initial = (
            context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoules_per_mole)
        )
        report["initial_energy_kj_mol"] = initial
        if not math.isfinite(initial):
            raise ValueError("Initial Amber19/GBn2 energy is nonfinite")
        LOGGER.info(
            "Starting Amber19SB/GBn2 minimization: %d movable atoms", len(movable)
        )
        mm.LocalEnergyMinimizer.minimize(
            context,
            tolerance=0.1 * unit.kilojoules_per_mole / unit.nanometer,
            maxIterations=20000,
        )
        state = context.getState(getPositions=True, getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        forces = state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoules_per_mole / unit.nanometer
        )
        if not math.isfinite(energy) or not all(
            math.isfinite(float(x)) for x in forces.flat
        ):
            raise ValueError("Final energy or forces nonfinite")
        report["final_energy_kj_mol"] = energy
        report["movable_rms_force_kj_mol_nm"] = float(
            math.sqrt(
                sum(float(x) ** 2 for i in movable for x in forces[i])
                / (3 * len(movable))
            )
        )
        report["force_tolerance_kj_mol_nm"] = 0.1
        report["max_iterations"] = 20000
        output = Modeller(modeller.topology, state.getPositions())
        output.delete([a for a in output.topology.atoms() if a.element.symbol == "H"])
        destination = output_dir / "refined_heavy.pdb"
        with destination.open("w") as handle:
            PDBFile.writeFile(output.topology, output.positions, handle, keepIds=True)
        report["output"] = str(destination)
        report["output_sha256"] = sha256(destination)
        report["geometry"] = geometry_audit(
            destination, original, dict(expected), termini
        )
        report["force_converged"] = report["movable_rms_force_kj_mol_nm"] <= 0.1
        report["state"] = (
            "accepted" if not report["geometry"]["failures"] else "rejected"
        )
        if report["state"] != "accepted":
            raise ValueError(
                f"Geometry acceptance failed: {report['geometry']['failures']}"
            )
        return report
    except Exception as exc:
        report["state"] = "rejected"
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        report_path.write_text(json.dumps(report, indent=2) + "\n")


def main() -> None:
    """Parse explicit candidate, original experimental core, and fresh output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--outputdir", type=Path, required=True)
    parser.add_argument("--termini", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    refine(
        args.candidate.resolve(),
        args.original.resolve(),
        args.outputdir.resolve(),
        args.termini.resolve(),
    )


if __name__ == "__main__":
    main()
