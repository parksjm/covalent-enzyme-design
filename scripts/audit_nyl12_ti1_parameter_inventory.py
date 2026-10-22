"""Inventory the current provisional TI1 parameters without accepting them."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

from parmed.amber import AmberParm
from parmed.topologyobjects import Atom


def main() -> None:
    """Retain explicit reactive atom/bond/angle/torsion data for validation.

    The input is the canonical 330-residue PA6 alpha/beta diagnostic topology.
    This read-only inventory never certifies a parameter set for refinement.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.output.exists():
        raise FileExistsError(args.output)
    model = AmberParm(str(args.topology))
    if len(model.residues) != 330:
        raise ValueError("Expected canonical alpha/beta diagnostic topology")
    thr, lig, water = [model.residues[i] for i in (232, 328, 329)]
    if [r.name for r in (thr, lig, water)] != ["THR", "LIG", "WAT"]:
        raise ValueError("Unexpected reactive residue identities")
    t, ligand = [{a.name: a for a in r.atoms} for r in (thr, lig)]
    og1, carbon = t["OG1"], ligand["C54"]
    if set(carbon.bond_partners) != {og1, ligand["O49"], ligand["N59"], ligand["C62"]}:
        raise ValueError("Topology does not have the intended TI1 connectivity")
    if {a.name for a in t["N"].bond_partners if a.atomic_number == 1} != {
        "H1",
        "H2",
        "H3",
    }:
        raise ValueError("Reacting N is not NH3+")
    selected = [*thr.atoms, *lig.atoms, *water.atoms]
    atoms = [
        {
            "index": a.idx + 1,
            "residue": a.residue.idx + 1,
            "name": a.name,
            "element_number": a.atomic_number,
            "type": a.type,
            "partial_charge": a.charge,
            "lj_rmin_half_A": a.rmin,
            "lj_epsilon_kcal_mol": a.epsilon,
        }
        for a in selected
    ]

    def label(atom: Atom) -> str:
        return f"{atom.residue.idx + 1}:{atom.name}"

    bonds = [
        {"atoms": [label(b.atom1), label(b.atom2)], "k": b.type.k, "req_A": b.type.req}
        for b in model.bonds
        if carbon in (b.atom1, b.atom2) or og1 in (b.atom1, b.atom2)
    ]
    angles = [
        {
            "atoms": [label(a.atom1), label(a.atom2), label(a.atom3)],
            "k": a.type.k,
            "equilibrium_degrees": a.type.theteq,
        }
        for a in model.angles
        if a.atom2 in (carbon, og1)
    ]
    torsions = []
    for term in model.dihedrals:
        if {term.atom2, term.atom3} != {og1, carbon}:
            continue
        types = list(term.type) if isinstance(term.type, list) else [term.type]
        torsions.append(
            {
                "atoms": [
                    label(a) for a in (term.atom1, term.atom2, term.atom3, term.atom4)
                ],
                "improper": term.improper,
                "terms": [
                    {"k": v.phi_k, "periodicity": v.per, "phase_degrees": v.phase}
                    for v in types
                ],
            }
        )
    report = {
        "topology": str(args.topology.absolute()),
        "topology_sha256": hashlib.sha256(args.topology.read_bytes()).hexdigest(),
        "audit_type": "read-only provisional parameter inventory",
        "classical_ti1_ready": False,
        "validated_parameter_set": False,
        "formal_state_declared": "reacting Thr NH3+, tetrahedral substrate oxyanion, neutral water; partial charges do not establish formal charges",
        "reactive_atoms": atoms,
        "bonds_at_reactive_center": bonds,
        "angles_at_reactive_center": angles,
        "torsions_across_junction": torsions,
        "missing_evidence": [
            "QM reference for the coupled TI1 adduct and its stated environment",
            "coupled-state constrained RESP fit and cap-transfer validation",
            "TI1-specific bonded parameters and independent geometry/curvature/torsion checks",
            "acceptance tolerances and final validation tied to exact topology digest",
        ],
        "parameters_modified": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    logging.info(
        "Inventoried %d reactive atoms, %d bonds, %d angles, %d junction torsions; TI1 validated=false",
        len(atoms),
        len(bonds),
        len(angles),
        len(torsions),
    )


if __name__ == "__main__":
    main()
