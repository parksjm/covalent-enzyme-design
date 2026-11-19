"""Build the explicitly bonded TI1 chemical reference graph."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from nyl12_covalent.chemistry import load_pa6
from nyl12_covalent.terminal_template import _build_template

LOGGER = logging.getLogger(__name__)
AttackCIP = Literal["R", "S"]


def build_ti1_reference_graph(
    pa6_mol2_path: Path, attacked_carbon_cip: AttackCIP
) -> tuple[Chem.Mol, dict[str, object]]:
    """Build and validate the 103-atom covalent TI1 graph.

    The returned molecule has no conformer. Its atom order is stable: neutral
    Thr-NME template atoms, source PA6 atoms, then explicit water O/H1/H2.
    The caller must choose the attacked-carbon CIP configuration explicitly.

    Args:
        pa6_mol2_path: Full, explicitly protonated zwitterionic PA6 L4 MOL2.
        attacked_carbon_cip: Required R or S configuration at PA6 C54.

    Returns:
        An RDKit molecule and metadata containing stable atom-origin/index
        mappings and hydrogen provenance.

    Raises:
        ValueError: If the source chemistry, requested stereochemistry, or
            resulting graph fails any explicit inventory or valence check.
    """
    if attacked_carbon_cip not in {"R", "S"}:
        raise ValueError("attacked_carbon_cip must be explicitly 'R' or 'S'")
    if not pa6_mol2_path.is_file():
        raise FileNotFoundError(f"PA6 MOL2 source does not exist: {pa6_mol2_path}")

    pa6 = load_pa6(pa6_mol2_path)
    validate_pa6_zwitterion(pa6.mol)
    template, _ = _build_template()
    combined = Chem.CombineMols(template, pa6.mol)
    rw = Chem.RWMol(combined)
    offset = template.GetNumAtoms()

    def index(residue: str, atom_name: str) -> int:
        if residue == "THR":
            matches = [
                atom.GetIdx()
                for atom in template.GetAtoms()
                if atom.GetProp("residue_name") == residue
                and atom.GetProp("atom_name") == atom_name
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected one template atom {residue}:{atom_name}; "
                    f"found {len(matches)}"
                )
            return matches[0]
        matches = [
            atom.GetIdx()
            for atom in pa6.mol.GetAtoms()
            if atom.GetProp("name") == atom_name
        ]
        if len(matches) != 1:
            raise ValueError(f"Expected one PA6 atom {atom_name}; found {len(matches)}")
        return offset + matches[0]

    og1 = index("THR", "OG1")
    thr_n = index("THR", "N")
    hg1 = index("THR", "HG1")
    c54 = index("PA6", "C54")
    o49 = index("PA6", "O49")
    n59 = index("PA6", "N59")
    c62 = index("PA6", "C62")

    _require_bond(rw, og1, hg1, rdchem.BondType.SINGLE, "Thr OG1-HG1")
    _require_bond(rw, c54, o49, rdchem.BondType.DOUBLE, "PA6 C54=O49")
    _require_bond(rw, c54, n59, rdchem.BondType.SINGLE, "PA6 C54-N59")
    _require_bond(rw, c54, c62, rdchem.BondType.SINGLE, "PA6 C54-C62")
    if rw.GetBondBetweenAtoms(thr_n, hg1) is not None:
        raise ValueError("Thr N already has a bond to HG1")

    rw.RemoveBond(og1, hg1)
    rw.AddBond(thr_n, hg1, rdchem.BondType.SINGLE)
    rw.GetAtomWithIdx(thr_n).SetFormalCharge(1)
    rw.GetAtomWithIdx(og1).SetFormalCharge(0)
    rw.GetAtomWithIdx(hg1).SetProp("atom_name", "H3")
    rw.GetAtomWithIdx(hg1).SetProp("_Name", "H3")
    rw.RemoveBond(c54, o49)
    rw.AddBond(c54, o49, rdchem.BondType.SINGLE)
    rw.GetAtomWithIdx(o49).SetFormalCharge(-1)
    rw.AddBond(og1, c54, rdchem.BondType.SINGLE)

    water_o = _add_water(rw)
    molecule = rw.GetMol()
    Chem.SanitizeMol(molecule)
    _set_cip(molecule, c54, attacked_carbon_cip)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    molecule.RemoveAllConformers()
    _annotate_origins(molecule, template.GetNumAtoms(), pa6.mol.GetNumAtoms())
    _validate_product(molecule, c54, attacked_carbon_cip, water_o)

    origins = {atom.GetProp("origin_id"): atom.GetIdx() for atom in molecule.GetAtoms()}
    atom_records = [
        {
            "index": atom.GetIdx(),
            "origin_id": atom.GetProp("origin_id"),
            "residue_id": atom.GetProp("residue_id"),
            "atom_name": atom.GetProp("atom_name"),
            "element": atom.GetSymbol(),
            "formal_charge": atom.GetFormalCharge(),
        }
        for atom in molecule.GetAtoms()
    ]
    metadata: dict[str, object] = {
        "formula": CalcMolFormula(molecule),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "attacked_carbon_cip": attacked_carbon_cip,
        "atom_index_by_origin": origins,
        "atoms": atom_records,
        "hydrogen_bookkeeping": {
            "origin_id": "THR:HG1",
            "source_atom_name": "HG1",
            "atom_name_after_transfer": "H3",
            "bonded_to_after_transfer": "THR:N",
            "identity_note": "Origin and H3 name document graph bookkeeping, not literal relay-proton identity.",
        },
    }
    LOGGER.info(
        "Built TI1 reference graph: %s, charge %d",
        metadata["formula"],
        metadata["formal_charge"],
    )
    return molecule, metadata


def validate_pa6_zwitterion(molecule: Chem.Mol) -> None:
    """Require the canonical PA6 L4 atom inventory and zwitterionic termini.

    Args:
        molecule: Explicit-hydrogen, named PA6 L4 molecular graph.

    Raises:
        ValueError: If atom count, net charge, required names, or terminal
            protonation differs from the canonical soluble-oligomer reference.
    """
    if molecule.GetNumAtoms() != 79 or Chem.GetFormalCharge(molecule) != 0:
        raise ValueError("PA6 source must be the 79-atom neutral zwitterion")
    named = {
        atom.GetProp("name"): atom
        for atom in molecule.GetAtoms()
        if atom.HasProp("name")
    }
    required = {"N56", "C53", "O48", "O51", "C54", "O49", "N59", "C62", "H3"}
    missing = required - named.keys()
    if missing:
        raise ValueError(f"PA6 source is missing named atoms: {sorted(missing)}")
    nitrogen = named["N56"]
    attached_h = [
        neighbor for neighbor in nitrogen.GetNeighbors() if neighbor.GetSymbol() == "H"
    ]
    if nitrogen.GetFormalCharge() != 1 or len(attached_h) != 3:
        raise ValueError("PA6 N56 must remain explicit NH3+")
    if named["H3"].GetIdx() not in {atom.GetIdx() for atom in attached_h}:
        raise ValueError("PA6 H3 provenance must remain attached to N56")
    acid_oxygens = [named[name] for name in ("O48", "O51")]
    if sum(atom.GetFormalCharge() for atom in acid_oxygens) != -1:
        raise ValueError("PA6 C53 terminus must be carboxylate")


def _require_bond(
    molecule: Chem.RWMol,
    first: int,
    second: int,
    bond_type: rdchem.BondType,
    label: str,
) -> None:
    """Require a specific existing bond before applying chemistry edits."""
    bond = molecule.GetBondBetweenAtoms(first, second)
    if bond is None or bond.GetBondType() != bond_type:
        raise ValueError(f"Expected {label} bond of type {bond_type}")


def _add_water(molecule: Chem.RWMol) -> int:
    """Append explicit neutral water atoms and return the oxygen index."""
    oxygen = Chem.Atom(8)
    oxygen_idx = molecule.AddAtom(oxygen)
    hydrogen_indices = []
    for _ in range(2):
        hydrogen_indices.append(molecule.AddAtom(Chem.Atom(1)))
    for hydrogen_idx in hydrogen_indices:
        molecule.AddBond(oxygen_idx, hydrogen_idx, rdchem.BondType.SINGLE)
    return oxygen_idx


def _set_cip(molecule: Chem.Mol, atom_index: int, target: AttackCIP) -> None:
    """Set a tetrahedral tag to the explicitly requested CIP label."""
    atom = molecule.GetAtomWithIdx(atom_index)
    atom.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    if not atom.HasProp("_CIPCode"):
        raise ValueError("Attacked carbon C54 is not stereogenic in the graph")
    if atom.GetProp("_CIPCode") != target:
        atom.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != target:
        raise ValueError(f"Could not encode requested C54 stereochemistry {target}")


def _annotate_origins(molecule: Chem.Mol, template_count: int, pa6_count: int) -> None:
    """Attach stable residue-qualified names and source indices to each atom."""
    for atom in molecule.GetAtoms():
        idx = atom.GetIdx()
        if idx < template_count:
            residue = atom.GetProp("residue_name")
            atom_name = atom.GetProp("atom_name")
            origin_name = (
                "HG1"
                if idx == _template_hg1_index(molecule, template_count)
                else atom_name
            )
            source_index = idx
        elif idx < template_count + pa6_count:
            source_index = idx - template_count
            atom_name = atom.GetProp("name")
            residue = "PA6"
            origin_name = atom_name
            atom.SetProp("atom_name", atom_name)
        else:
            source_index = -1
            residue = "WAT"
            atom_name = (
                "O"
                if atom.GetSymbol() == "O"
                else f"H{idx - template_count - pa6_count}"
            )
            origin_name = atom_name
            atom.SetProp("atom_name", atom_name)
        atom.SetProp("residue_id", residue)
        atom.SetProp("origin_id", f"{residue}:{origin_name}")
        atom.SetIntProp("source_atom_index", source_index)


def _template_hg1_index(molecule: Chem.Mol, template_count: int) -> int:
    """Find the transferred template proton by its stable origin property."""
    matches = [
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetIdx() < template_count
        and atom.HasProp("residue_name")
        and atom.GetProp("residue_name") == "THR"
        and atom.HasProp("atom_name")
        and atom.GetProp("atom_name") == "H3"
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one transferred Thr H3; found {len(matches)}")
    return matches[0]


def _validate_product(
    molecule: Chem.Mol, c54: int, target: AttackCIP, water_o: int
) -> None:
    """Check final atom, charge, formula, connectivity, valence and stereo."""
    if molecule.GetNumAtoms() != 103:
        raise ValueError(f"TI1 graph must have 103 atoms, got {molecule.GetNumAtoms()}")
    if CalcMolFormula(molecule) != "C29H60N6O8":
        raise ValueError(f"Unexpected TI1 formula: {CalcMolFormula(molecule)}")
    if Chem.GetFormalCharge(molecule) != 0:
        raise ValueError("TI1 reference graph must be net neutral")
    named = {atom.GetProp("origin_id"): atom for atom in molecule.GetAtoms()}
    expected_bonds = (
        ("THR:OG1", "PA6:C54", rdchem.BondType.SINGLE),
        ("PA6:C54", "PA6:O49", rdchem.BondType.SINGLE),
        ("PA6:C54", "PA6:N59", rdchem.BondType.SINGLE),
        ("THR:N", "THR:HG1", rdchem.BondType.SINGLE),
        ("WAT:O", "WAT:H1", rdchem.BondType.SINGLE),
        ("WAT:O", "WAT:H2", rdchem.BondType.SINGLE),
    )
    for first, second, bond_type in expected_bonds:
        bond = molecule.GetBondBetweenAtoms(
            named[first].GetIdx(), named[second].GetIdx()
        )
        if bond is None or bond.GetBondType() != bond_type:
            raise ValueError(f"Invalid TI1 connectivity: {first}-{second}")
    if molecule.GetBondBetweenAtoms(
        named["THR:OG1"].GetIdx(), named["THR:HG1"].GetIdx()
    ):
        raise ValueError("Transferred Thr H3 remains bonded to OG1")
    if molecule.GetAtomWithIdx(c54).GetProp("_CIPCode") != target:
        raise ValueError("C54 CIP assignment differs from the requested choice")
    if named["THR:OG1"].GetFormalCharge() != 0 or named["THR:N"].GetFormalCharge() != 1:
        raise ValueError("Thr protonation state is incorrect")
    if named["PA6:O49"].GetFormalCharge() != -1:
        raise ValueError("PA6 O49 must carry the oxyanion charge")
    if molecule.GetAtomWithIdx(water_o).GetFormalCharge() != 0:
        raise ValueError("Added water must be neutral")
    if molecule.GetNumConformers():
        raise ValueError(
            "Reference graph must not contain generated or stale coordinates"
        )
