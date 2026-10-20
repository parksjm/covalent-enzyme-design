"""Repair a specific LEaP exclusion-list defect and audit the acyl bond."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from parmed.amber import AmberFormat, AmberParm


def validate(output: Path, carbonyl_name: str) -> dict[str, object]:
    """Validate the generated Amber topology and write a machine-readable audit.

    AmberTools26 LEaP inserts one out-of-range exclusion integer at the new
    cross-residue bond in this system. The repair removes exactly that extra
    integer; any other pattern fails instead of changing topology silently.

    Args:
        output: Acyl model directory containing LEaP outputs.
        carbonyl_name: Stable ligand atom name of the ester carbonyl carbon.

    Returns:
        Atom, bond, charge, and exclusion-repair measurements.
    """
    raw = AmberFormat(str(output / "acyl.raw.prmtop"))
    atom_count = raw.parm_data["POINTERS"][0]
    exclusions = raw.parm_data["EXCLUDED_ATOMS_LIST"]
    expected_exclusions = sum(raw.parm_data["NUMBER_EXCLUDED_ATOMS"])
    bad = [index for index, atom in enumerate(exclusions) if atom < 0 or atom > atom_count]
    if bad:
        if len(bad) != 1 or len(exclusions) != expected_exclusions + 1:
            raise ValueError(f"Unexpected LEaP exclusion corruption: {bad}")
        exclusions.pop(bad[0])
    if len(exclusions) != expected_exclusions:
        raise ValueError("Amber exclusion count does not match its pointer table")
    raw.parm_data["EXCLUDED_ATOMS_LIST"] = exclusions
    clean_path = output / "acyl.prmtop"
    raw.write_parm(str(clean_path))
    topology = AmberParm(str(clean_path))
    if len(topology.atoms) != atom_count or len(topology.residues) != 1313:
        raise ValueError("Unexpected mature Nyl12 atom or residue inventory")
    thr = topology.residues[928]
    ligand = topology.residues[1312]
    if thr.name != "THR" or ligand.name != "LIG":
        raise ValueError("Thr233/LIG residue order changed")
    thr_atoms = {atom.name: atom for atom in thr.atoms}
    lig_atoms = {atom.name: atom for atom in ligand.atoms}
    if "H3" in thr_atoms or "HG1" in thr_atoms:
        raise ValueError("Catalytic Thr still bears transferred hydrogens")
    if carbonyl_name not in lig_atoms:
        raise ValueError(f"Missing selected carbonyl atom {carbonyl_name}")
    oxygen = thr_atoms["OG1"]
    carbon = lig_atoms[carbonyl_name]
    if carbon not in oxygen.bond_partners:
        raise ValueError("Thr OG1–acyl carbonyl bond is absent")
    if len(carbon.bond_partners) != 3:
        raise ValueError("Acyl carbonyl carbon lacks three bonded neighbors")
    total_charge = sum(atom.charge for atom in topology.atoms)
    qm_charge = sum(atom.charge for atom in [*thr.atoms, *ligand.atoms])
    audit: dict[str, object] = {
        "atom_count": atom_count,
        "residue_count": len(topology.residues),
        "retained_ligand_atoms": len(ligand.atoms),
        "total_charge": total_charge,
        "qm_region_partial_charge": qm_charge,
        "qm_region_formal_charge": 1,
        "carbonyl_partial_charge": carbon.charge,
        "ester_bond": ["929:OG1", f"1313:{carbonyl_name}"],
        "leap_exclusion_extra_entries_removed": len(bad),
        "leap_exclusion_entry_index": bad[0] if bad else None,
    }
    (output / "topology_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit


if __name__ == "__main__":
    validate(Path(sys.argv[1]), sys.argv[2])
