"""Audit the tetrahedral adduct, cationic beta N-terminus, and relay water."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from parmed.amber import AmberFormat, AmberParm


def validate(
    output: Path, carbonyl_name: str, oxygen_name: str, leaving_n: str
) -> dict[str, object]:
    """Require the first covalent intermediate's exact atom and bond inventory.

    Args:
        output: TI1 LEaP build directory.
        carbonyl_name: Attacked substrate carbonyl carbon.
        oxygen_name: Former carbonyl oxygen, now an oxyanion.
        leaving_n: Amide nitrogen that must remain bonded in TI1.

    Returns:
        Topology and protonation audit saved as JSON.
    """
    raw = AmberFormat(str(output / "ti1.raw.prmtop"))
    atom_count = raw.parm_data["POINTERS"][0]
    exclusions = raw.parm_data["EXCLUDED_ATOMS_LIST"]
    expected_exclusions = sum(raw.parm_data["NUMBER_EXCLUDED_ATOMS"])
    bad = [index for index, atom in enumerate(exclusions) if atom < 0 or atom > atom_count]
    if bad:
        if len(bad) != 1 or len(exclusions) != expected_exclusions + 1:
            raise ValueError(f"Unexpected LEaP exclusion corruption: {bad}")
        exclusions.pop(bad[0])
    if len(exclusions) != expected_exclusions:
        raise ValueError("TI1 exclusion count and pointer table disagree")
    raw.parm_data["EXCLUDED_ATOMS_LIST"] = exclusions
    topology_path = output / "ti1.prmtop"
    raw.write_parm(str(topology_path))
    topology = AmberParm(str(topology_path))
    if len(topology.residues) != 1314:
        raise ValueError("Expected mature Nyl12, intact ligand, and one relay water")
    thr, ligand, water = (topology.residues[index] for index in (928, 1312, 1313))
    if (thr.name, ligand.name, water.name) != ("THR", "LIG", "WAT"):
        raise ValueError("Catalytic Thr, substrate, or relay-water residue changed")
    thr_atoms = {atom.name: atom for atom in thr.atoms}
    lig_atoms = {atom.name: atom for atom in ligand.atoms}
    water_atoms = {atom.name: atom for atom in water.atoms}
    if set(water_atoms) != {"O", "H1", "H2"}:
        raise ValueError("Relay water must remain neutral TIP3P H2O")
    if "HG1" in thr_atoms or not {"H1", "H2", "H3"}.issubset(thr_atoms):
        raise ValueError("TI1 needs deprotonated OG1 and cationic beta N-H3")
    beta_n = thr_atoms["N"]
    if not all(thr_atoms[name] in beta_n.bond_partners for name in ("H1", "H2", "H3")):
        raise ValueError("Beta N-terminus does not bear all three hydrogens")
    for index in (1024, 1120, 1216):
        other = topology.residues[index]
        atoms = {atom.name: atom for atom in other.atoms}
        if other.name != "THR" or "H3" in atoms:
            raise ValueError("Nonreactive beta N terminus is not neutral NH2")
        if not all(atoms[name] in atoms["N"].bond_partners for name in ("H1", "H2")):
            raise ValueError("Nonreactive beta N terminus lacks NH2 bonds")
        if abs(sum(atom.charge for atom in other.atoms)) > 0.002:
            raise ValueError("Nonreactive beta N terminus has nonzero residue charge")
    carbon = lig_atoms[carbonyl_name]
    oxygen = lig_atoms[oxygen_name]
    nitrogen = lig_atoms[leaving_n]
    if len(carbon.bond_partners) != 4:
        raise ValueError("TI1 carbonyl carbon is not tetra-coordinate")
    if not all(
        atom in carbon.bond_partners
        for atom in (thr_atoms["OG1"], oxygen, nitrogen)
    ):
        raise ValueError("TI1 lost ester, oxyanion, or intact amide C-N bond")
    if len(oxygen.bond_partners) != 1:
        raise ValueError("Oxyanion should be singly bonded to carbonyl carbon")
    if abs(sum(atom.charge for atom in water.atoms)) > 0.001:
        raise ValueError("Relay water is not neutral")
    qm_atoms = [*thr.atoms, *ligand.atoms, *water.atoms]
    audit: dict[str, object] = {
        "atom_count": len(topology.atoms),
        "residue_count": len(topology.residues),
        "ligand_atoms": len(ligand.atoms),
        "total_charge": sum(atom.charge for atom in topology.atoms),
        "qm_region_partial_charge": sum(atom.charge for atom in qm_atoms),
        "qm_region_formal_charge": 0,
        "beta_n_terminus": "cationic NH3+ (H1,H2,H3 bonded)",
        "other_beta_n_termini": "neutral NH2 (H1,H2 bonded; H3 absent)",
        "og1_hg1": "absent, proton accounted for at beta N H3",
        "relay_water": "neutral H2O",
        "oxyanion_partial_charge": oxygen.charge,
        "tetrahedral_carbon_bonded_atoms": [atom.name for atom in carbon.bond_partners],
        "leap_exclusion_extra_entries_removed": len(bad),
    }
    (output / "ti1_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit


if __name__ == "__main__":
    validate(Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4])
