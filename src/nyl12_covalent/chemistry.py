"""Build assayed Nyl12 substrate graphs and verify hydrolysis products."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ASSAY_PRODUCTS = {
    "pa6_l2": ("L2", "C12H24N2O3"),
    "poly_xama_maxa": ("MAXA", "C26H42N4O5"),
}

MONADS = {
    "M": "[NH2:1]CCCCCC[NH2:2]",
    "X": "[NH2:1]Cc1ccc(C[NH2:2])cc1",
    "A": "[OH:1]C(=O)CCCCC(=O)[OH:2]",
}


@dataclass(frozen=True)
class Substrate:
    """A substrate with a named, product-defining amide bond.

    Attributes:
        key: Stable substrate target identifier.
        mol: Explicit-hydrogen RDKit graph with named atoms.
        carbonyl_c: Name of the attacked carbonyl carbon.
        carbonyl_o: Name of its oxygen.
        leaving_n: Name of the scissile amide nitrogen.
        source: Chemistry provenance.
    """

    key: str
    mol: Chem.Mol
    carbonyl_c: str
    carbonyl_o: str
    leaving_n: str
    source: str


def atom_by_name(mol: Chem.Mol, name: str) -> Chem.Atom:
    """Return the unique atom bearing a stable name.

    Args:
        mol: Named molecular graph.
        name: Required atom name.

    Raises:
        ValueError: If the name is absent or ambiguous.
    """
    matches = [a for a in mol.GetAtoms() if a.HasProp("name") and a.GetProp("name") == name]
    if len(matches) != 1:
        raise ValueError(f"Expected one atom named {name}; found {len(matches)}")
    return matches[0]


def load_pa6(mol2_path: Path) -> Substrate:
    """Load the original PA6 L4z model with its affinity-mode atom names.

    Args:
        mol2_path: GAFF input MOL2 associated with the selected WT model.

    Returns:
        Named PA6 substrate with the L2-producing amide selected.
    """
    mol = Chem.MolFromMol2File(str(mol2_path), removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError(f"Cannot parse source PA6 MOL2: {mol2_path}")
    for atom in mol.GetAtoms():
        if not atom.HasProp("_TriposAtomName"):
            raise ValueError("Source PA6 MOL2 lost an atom name")
        atom.SetProp("name", atom.GetProp("_TriposAtomName"))
    substrate = Substrate("pa6_l2", mol, "C54", "O49", "N59", str(mol2_path))
    validate_amide(substrate)
    return substrate


def _tagged_monad(kind: str, position: int) -> Chem.Mol:
    """Construct one repeat monad with tagged functional-group ends."""
    mol = Chem.MolFromSmiles(MONADS[kind])
    if mol is None:
        raise ValueError(f"Invalid {kind} monad template")
    for atom in mol.GetAtoms():
        atom.SetIntProp("monad_index", position)
        atom.SetProp("monad_kind", kind)
    return mol


def _tagged_atom(mol: Chem.Mol, position: int, map_number: int) -> Chem.Atom:
    """Select a monad endpoint by position and atom-map number."""
    matches = [
        atom
        for atom in mol.GetAtoms()
        if atom.HasProp("monad_index")
        and atom.GetIntProp("monad_index") == position
        and atom.GetAtomMapNum() == map_number
    ]
    if len(matches) != 1:
        raise ValueError(f"Ambiguous endpoint {position}:{map_number}")
    return matches[0]


def _couple(left: Chem.Mol, right: Chem.Mol, position: int, left_kind: str) -> Chem.Mol:
    """Form one amide, removing the relevant acid hydroxyl oxygen."""
    combined = Chem.CombineMols(left, right)
    acid_position = position - 1 if left_kind == "A" else position
    acid_end = 2 if left_kind == "A" else 1
    amine_position = position if left_kind == "A" else position - 1
    amine_end = 1 if left_kind == "A" else 2
    acid_oh = _tagged_atom(combined, acid_position, acid_end)
    amine = _tagged_atom(combined, amine_position, amine_end)
    if acid_oh.GetAtomicNum() != 8 or amine.GetAtomicNum() != 7:
        raise ValueError("Monads do not present an acid/amine junction")
    neighbors = list(acid_oh.GetNeighbors())
    if len(neighbors) != 1 or neighbors[0].GetAtomicNum() != 6:
        raise ValueError("Acid endpoint lacks a unique carbonyl carbon")
    carbon = neighbors[0]
    carbon.SetProp("junction_to", f"{position - 1}|{position}")
    amine.SetProp("junction_to", f"{position - 1}|{position}")
    rw = Chem.RWMol(combined)
    rw.RemoveAtom(acid_oh.GetIdx())
    carbon_index = carbon.GetIdx() - int(acid_oh.GetIdx() < carbon.GetIdx())
    amine_index = amine.GetIdx() - int(acid_oh.GetIdx() < amine.GetIdx())
    rw.AddBond(carbon_index, amine_index, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(amine_index).SetNumExplicitHs(1)
    result = rw.GetMol()
    Chem.SanitizeMol(result)
    return result


def build_poly_xama_maxa() -> Substrate:
    """Build poly(X,AMA)'s MAXA-producing A|M dimer register.

    Returns:
        Zwitterionic M-A-X-A-M-A-X-A oligomer with a central A|M cut.
    """
    sequence = "MAXAMAXA"
    mol = _tagged_monad(sequence[0], 1)
    for position, kind in enumerate(sequence[1:], start=2):
        mol = _couple(mol, _tagged_monad(kind, position), position, sequence[position - 2])
    first_n = _tagged_atom(mol, 1, 1)
    last_o = _tagged_atom(mol, len(sequence), 2)
    first_n.SetFormalCharge(1)
    first_n.SetNumExplicitHs(3)
    last_o.SetFormalCharge(-1)
    last_o.SetNumExplicitHs(0)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    for atom in mol.GetAtoms():
        atom.SetProp("name", f"{atom.GetSymbol().upper()}{atom.GetIdx() + 1:03d}")
    junction = [
        atom for atom in mol.GetAtoms() if atom.HasProp("junction_to") and atom.GetProp("junction_to") == "4|5"
    ]
    carbon = next((a for a in junction if a.GetAtomicNum() == 6), None)
    nitrogen = next((a for a in junction if a.GetAtomicNum() == 7), None)
    if carbon is None or nitrogen is None:
        raise ValueError("Central A|M junction is incomplete")
    oxygen = next(
        (
            bond.GetOtherAtom(carbon)
            for bond in carbon.GetBonds()
            if bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(carbon).GetAtomicNum() == 8
        ),
        None,
    )
    if oxygen is None:
        raise ValueError("Central A|M junction lacks a carbonyl oxygen")
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        if not atom.HasProp("name"):
            atom.SetProp("name", f"H{atom.GetIdx() + 1:03d}")
    substrate = Substrate(
        "poly_xama_maxa", mol, carbon.GetProp("name"), oxygen.GetProp("name"),
        nitrogen.GetProp("name"), "Dishner 2026 JACS; M-A-X-A repeated twice",
    )
    validate_amide(substrate)
    return substrate


def validate_amide(substrate: Substrate) -> None:
    """Reject a scissile selection inconsistent with an intact amide."""
    carbon = atom_by_name(substrate.mol, substrate.carbonyl_c)
    oxygen = atom_by_name(substrate.mol, substrate.carbonyl_o)
    nitrogen = atom_by_name(substrate.mol, substrate.leaving_n)
    co = substrate.mol.GetBondBetweenAtoms(carbon.GetIdx(), oxygen.GetIdx())
    cn = substrate.mol.GetBondBetweenAtoms(carbon.GetIdx(), nitrogen.GetIdx())
    if co is None or co.GetBondType() != Chem.BondType.DOUBLE:
        raise ValueError("Selected carbonyl C/O is not a double bond")
    if cn is None or cn.GetBondType() != Chem.BondType.SINGLE:
        raise ValueError("Selected C/N is not an amide single bond")


def _hydrolysis_products(substrate: Substrate) -> tuple[Chem.Mol, ...]:
    """Hydrolyze the selected amide with explicit water-atom bookkeeping."""
    mol = Chem.RemoveHs(substrate.mol)
    carbon = atom_by_name(mol, substrate.carbonyl_c)
    nitrogen = atom_by_name(mol, substrate.leaving_n)
    rw = Chem.RWMol(mol)
    rw.RemoveBond(carbon.GetIdx(), nitrogen.GetIdx())
    oxygen = Chem.Atom(8)
    oxygen.SetFormalCharge(-1)
    oxygen.SetProp("name", "OWAT")
    oxygen_idx = rw.AddAtom(oxygen)
    rw.AddBond(carbon.GetIdx(), oxygen_idx, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(nitrogen.GetIdx()).SetFormalCharge(1)
    rw.GetAtomWithIdx(nitrogen.GetIdx()).SetNumExplicitHs(3)
    hydrolyzed = rw.GetMol()
    Chem.SanitizeMol(hydrolyzed)
    return Chem.GetMolFrags(hydrolyzed, asMols=True, sanitizeFrags=True)


def product_report(substrate: Substrate) -> dict[str, object]:
    """Return a strict product identity report for an assayed Nyl12 target.

    Raises:
        ValueError: If hydrolysis does not generate the required product.
    """
    expected_name, expected_formula = ASSAY_PRODUCTS[substrate.key]
    products = _hydrolysis_products(substrate)
    formulas = [rdMolDescriptors.CalcMolFormula(product) for product in products]
    if len(products) != 2 or formulas.count(expected_formula) < 1:
        raise ValueError(
            f"{substrate.key} cleavage yields {formulas}; expected {expected_name} "
            f"({expected_formula})"
        )
    if substrate.key == "pa6_l2" and formulas.count(expected_formula) != 2:
        raise ValueError(f"PA6 L4 cleavage does not make two L2 fragments: {formulas}")
    return {
        "substrate": substrate.key,
        "assayed_product": expected_name,
        "expected_formula": expected_formula,
        "fragment_formulas": formulas,
        "fragment_smiles": [Chem.MolToSmiles(product) for product in products],
        "product_exact_mass": Descriptors.ExactMolWt(products[formulas.index(expected_formula)]),
        "scissile": {
            "carbonyl_c": substrate.carbonyl_c,
            "carbonyl_o": substrate.carbonyl_o,
            "leaving_n": substrate.leaving_n,
        },
        "source": substrate.source,
    }


def retained_names(substrate: Substrate) -> frozenset[str]:
    """Return all atoms on the carbonyl side of the selected amide cut."""
    carbon = atom_by_name(substrate.mol, substrate.carbonyl_c)
    nitrogen = atom_by_name(substrate.mol, substrate.leaving_n)
    blocked = frozenset((carbon.GetIdx(), nitrogen.GetIdx()))
    seen = {carbon.GetIdx()}
    pending = [carbon]
    while pending:
        for neighbor in pending.pop().GetNeighbors():
            if frozenset((carbon.GetIdx(), neighbor.GetIdx())) == blocked:
                continue
            if neighbor.GetIdx() not in seen:
                seen.add(neighbor.GetIdx())
                pending.append(neighbor)
    return frozenset(substrate.mol.GetAtomWithIdx(index).GetProp("name") for index in seen)
