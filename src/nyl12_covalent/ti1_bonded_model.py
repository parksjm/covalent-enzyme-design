"""Build a capped, local TI1 surrogate for bonded-parameter fitting."""

from __future__ import annotations

import logging
from typing import Literal

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from nyl12_covalent.terminal_template import _build_template

LOGGER = logging.getLogger(__name__)
AttackCIP = Literal["R", "S"]
_EXPECTED_THR_CIP = {"CA": "S", "CB": "R"}


def build_ti1_bonded_model(
    attacked_carbon_cip: AttackCIP,
) -> tuple[Chem.Mol, dict[str, object]]:
    """Build the ACE-capped local TI1 bonded-fitting surrogate.

    The molecule represents ACE-THR(-O-C*(O-)(CH3)(NHCH3))-NME. The
    transferred Thr hydroxyl proton is omitted into an external reservoir.
    This charge -1 truncated model is only for bonded fitting; it is not a
    native NH3+ model or a native ESP/RESP charge-transfer model.

    Args:
        attacked_carbon_cip: Required R or S configuration at TIL:C*.

    Returns:
        An RDKit molecule and metadata with stable atom-origin/index mappings.

    Raises:
        ValueError: If the requested stereochemistry or molecular invariants
            cannot be represented by the final graph.
    """
    if attacked_carbon_cip not in {"R", "S"}:
        raise ValueError("attacked_carbon_cip must be explicitly 'R' or 'S'")
    template, _ = _build_template()
    rw = Chem.RWMol(template)
    _prepare_template_atoms(rw)

    remove_indices = [
        _find_atom(rw, "THR", "H2"),
        _find_atom(rw, "THR", "HG1"),
    ]
    for atom_index in sorted(remove_indices, reverse=True):
        rw.RemoveAtom(atom_index)
    _rename_atom(rw, "THR", "H1", "H")

    ace_c = _add_atom(rw, "C", "ACE", "C")
    ace_o = _add_atom(rw, "O", "ACE", "O")
    ace_ch3 = _add_atom(rw, "C", "ACE", "CH3")
    til_c = _add_atom(rw, "C", "TIL", "C*")
    til_o = _add_atom(rw, "O", "TIL", "O", formal_charge=-1)
    til_n = _add_atom(rw, "N", "TIL", "N")
    til_acyl_ch3 = _add_atom(rw, "C", "TIL", "acylCH3")
    til_leaving_ch3 = _add_atom(rw, "C", "TIL", "leavingCH3")
    thr_n = _find_atom(rw, "THR", "N")
    thr_og1 = _find_atom(rw, "THR", "OG1")

    _add_bond(rw, ace_c, ace_o, rdchem.BondType.DOUBLE)
    _add_bond(rw, ace_c, ace_ch3, rdchem.BondType.SINGLE)
    _add_bond(rw, ace_c, thr_n, rdchem.BondType.SINGLE)
    _add_bond(rw, til_c, til_o, rdchem.BondType.SINGLE)
    _add_bond(rw, til_c, til_n, rdchem.BondType.SINGLE)
    _add_bond(rw, til_c, til_acyl_ch3, rdchem.BondType.SINGLE)
    _add_bond(rw, til_c, thr_og1, rdchem.BondType.SINGLE)
    _add_bond(rw, til_n, til_leaving_ch3, rdchem.BondType.SINGLE)

    for name, carbon_index in (
        ("ACE:CH3", ace_ch3),
        ("TIL:acylCH3", til_acyl_ch3),
        ("TIL:leavingCH3", til_leaving_ch3),
    ):
        _add_hydrogens(rw, carbon_index, name, 3)
    _add_hydrogen(rw, til_n, "TIL", "HN")

    molecule = rw.GetMol()
    Chem.SanitizeMol(molecule)
    _set_cip(molecule, til_c, attacked_carbon_cip)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    _validate_model(molecule, til_c, attacked_carbon_cip)

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
        "model_scope": "bonded_fit_surrogate_only",
        "native_esp_resp_transfer_eligible": False,
        "formula": CalcMolFormula(molecule),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "electron_count": _electron_count(molecule),
        "multiplicity": 1,
        "attacked_carbon_cip": attacked_carbon_cip,
        "atom_index_by_origin": origins,
        "atoms": atom_records,
        "proton_reservoir": {
            "omitted_origin_id": "THR:HG1",
            "note": "This truncation omits the hydroxyl proton into an external reservoir; it does not specify a relay identity.",
        },
    }
    LOGGER.info(
        "Built local TI1 bonded surrogate: %s, charge %d, %d electrons",
        metadata["formula"],
        metadata["formal_charge"],
        metadata["electron_count"],
    )
    return molecule, metadata


def _prepare_template_atoms(molecule: Chem.RWMol) -> None:
    """Add residue-qualified stable origin labels to template atoms."""
    for atom in molecule.GetAtoms():
        residue = atom.GetProp("residue_name")
        name = atom.GetProp("atom_name")
        atom.SetProp("residue_id", residue)
        atom.SetProp("origin_id", f"{residue}:{name}")


def _find_atom(molecule: Chem.RWMol, residue: str, atom_name: str) -> int:
    """Return the unique atom matching a template residue and atom name."""
    matches = [
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetProp("residue_name") == residue
        and atom.GetProp("atom_name") == atom_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one atom {residue}:{atom_name}; found {len(matches)}"
        )
    return matches[0]


def _rename_atom(
    molecule: Chem.RWMol, residue: str, old_name: str, new_name: str
) -> None:
    """Rename one template atom and update its stable origin label."""
    atom = molecule.GetAtomWithIdx(_find_atom(molecule, residue, old_name))
    atom.SetProp("atom_name", new_name)
    atom.SetProp("_Name", new_name)
    atom.SetProp("origin_id", f"{residue}:{new_name}")


def _add_atom(
    molecule: Chem.RWMol,
    symbol: str,
    residue: str,
    name: str,
    formal_charge: int = 0,
) -> int:
    """Append one named heavy atom and return its index."""
    atom = Chem.Atom(symbol)
    atom.SetFormalCharge(formal_charge)
    atom.SetProp("residue_name", residue)
    atom.SetProp("residue_id", residue)
    atom.SetProp("atom_name", name)
    atom.SetProp("_Name", name)
    atom.SetProp("origin_id", f"{residue}:{name}")
    return molecule.AddAtom(atom)


def _add_hydrogen(
    molecule: Chem.RWMol, parent_index: int, residue: str, name: str
) -> int:
    """Add one named explicit hydrogen and connect it to its parent."""
    index = _add_atom(molecule, "H", residue, name)
    molecule.AddBond(parent_index, index, rdchem.BondType.SINGLE)
    return index


def _add_hydrogens(
    molecule: Chem.RWMol, parent_index: int, prefix: str, count: int
) -> None:
    """Add a deterministic numbered set of explicit methyl hydrogens."""
    residue, label = prefix.split(":", maxsplit=1)
    for position in range(1, count + 1):
        _add_hydrogen(molecule, parent_index, residue, f"{label}H{position}")


def _add_bond(
    molecule: Chem.RWMol, first: int, second: int, bond_type: rdchem.BondType
) -> None:
    """Add a bond after rejecting duplicate connectivity."""
    if molecule.GetBondBetweenAtoms(first, second) is not None:
        raise ValueError(f"Duplicate bond between atoms {first} and {second}")
    molecule.AddBond(first, second, bond_type)


def _set_cip(molecule: Chem.Mol, atom_index: int, target: AttackCIP) -> None:
    """Set tetrahedral stereo to an explicit requested CIP configuration."""
    atom = molecule.GetAtomWithIdx(atom_index)
    atom.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    if not atom.HasProp("_CIPCode"):
        raise ValueError("TIL:C* is not stereogenic")
    if atom.GetProp("_CIPCode") != target:
        atom.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != target:
        raise ValueError(f"Could not encode requested TIL:C* stereochemistry {target}")


def _electron_count(molecule: Chem.Mol) -> int:
    """Return the total electron count implied by nuclei and formal charge."""
    return sum(
        atom.GetAtomicNum() for atom in molecule.GetAtoms()
    ) - Chem.GetFormalCharge(molecule)


def _validate_model(
    molecule: Chem.Mol, attacked_carbon: int, target: AttackCIP
) -> None:
    """Reject graphs inconsistent with the declared truncated-model scope."""
    if molecule.GetNumAtoms() != 37 or CalcMolFormula(molecule) != "C10H20N3O4-":
        raise ValueError("Bonded surrogate must contain C10H20N3O4- (37 atoms)")
    if Chem.GetFormalCharge(molecule) != -1 or _electron_count(molecule) != 134:
        raise ValueError("Bonded surrogate must be a 134-electron anion")
    if molecule.GetNumConformers() != 0:
        raise ValueError("Bonded surrogate graph must not contain coordinates")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    for name, cip in _EXPECTED_THR_CIP.items():
        atom = next(
            atom
            for atom in molecule.GetAtoms()
            if atom.GetProp("origin_id") == f"THR:{name}"
        )
        if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != cip:
            raise ValueError(f"THR:{name} stereochemistry must remain {cip}")
    carbon = molecule.GetAtomWithIdx(attacked_carbon)
    if not carbon.HasProp("_CIPCode") or carbon.GetProp("_CIPCode") != target:
        raise ValueError("Attacked-carbon stereochemistry does not match the request")
