"""Authentic noncovalent E-S ligand CCD and predictor-specific atom mappings."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

ROOT = Path("/data/jpc/covalent-enzyme-design")
CONFIG = ROOT / "design/nyl12_ligandmpnn_pilot/boltz_validation/configuration"
SMILES = "[NH3+]CCCCCC(NCCCCCC(NCCCCCC(NCCCCCC([O-])=O)=O)=O)=O"
COMPONENT = "L4ESB26"


def digest(path: Path) -> str:
    """Return the SHA256 identity of an existing file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_molecule(path: Path) -> Any:
    """Read and strictly validate the authentic79-atom zwitterionic L4 graph.

    Args:
        path: Existing validated whole-L4 Mol2 source.

    Returns:
        RDKit molecule retaining physical atom names and formal charges.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    mol = Chem.MolFromMol2File(str(path), removeHs=False)
    if mol is None or mol.GetNumAtoms() != 79 or Chem.GetFormalCharge(mol) != 0:
        raise ValueError("Expected79physicalatomsandnet-neutralzwitterionicL4")
    heavy = Chem.RemoveHs(mol)
    expected = Chem.MolFromSmiles(SMILES)
    if heavy.GetNumAtoms() != 33 or Chem.MolToSmiles(heavy) != Chem.MolToSmiles(
        expected
    ):
        raise ValueError("SourcegraphisnotvalidatednoncovalentzwitterionicL4")
    names = [a.GetProp("_TriposAtomName") for a in mol.GetAtoms()]
    if len(set(names)) != len(names):
        raise ValueError("Sourceatomnamesmustbeunique")
    return mol


def ligand_maps(mol: Any, native_table: dict) -> dict[str, dict[str, str]]:
    """Map actual predictor atom names onto original source-Mol2 names.

    Args:
        mol: Strictly validated source molecule with named atoms.
        native_table: Atom/bond table emitted by installed Boltz's actual parser.

    Returns:
        AF3 and Boltz dictionaries mapping model atom name to source Mol2 name.
    """
    if native_table["smiles"] != SMILES:
        raise ValueError("NativeBoltzatomtableusesdifferentSMILES")
    heavy = Chem.RemoveHs(mol)
    query = Chem.MolFromSmiles(SMILES)
    matches = heavy.GetSubstructMatches(query, uniquify=False)
    if len(matches) != 1:
        raise ValueError("SMILEStoMol2heavyatommapmustbeunambiguous")
    match = matches[0]
    atoms = native_table["atoms"]
    if len(atoms) != 33 or {a["index0based"] for a in atoms} != set(range(33)):
        raise ValueError("NativeBoltztablemustcontain33uniqueheavyatomindices")
    names = {}
    for atom in atoms:
        index = atom["index0based"]
        q = query.GetAtomWithIdx(index)
        if (atom["element"], atom["formal_charge"]) != (
            q.GetSymbol(),
            q.GetFormalCharge(),
        ):
            raise ValueError("NativeBoltzSMILESatomorderorformalchargesdiffer")
        names[atom["name"]] = heavy.GetAtomWithIdx(match[index]).GetProp(
            "_TriposAtomName"
        )
    edges = {
        tuple(sorted((b["first0based"], b["second0based"]))): b["order"]
        for b in native_table["bonds"]
    }
    expected = {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))): b.GetBondTypeAsDouble()
        for b in query.GetBonds()
    }
    if edges != expected or len(names) != 33 or len(set(names.values())) != 33:
        raise ValueError("NativeBoltzbondgraphoratomnamemapchanged")
    identity = {
        a.GetProp("_TriposAtomName"): a.GetProp("_TriposAtomName")
        for a in heavy.GetAtoms()
    }
    return dict(af3=identity, boltz=names)


def write_es_ccd(mol: Any, path: Path) -> None:
    """Write a fresh ideal-conformer E-S CCD without protein coordinates or bonds.

    Args:
        mol: Authentic full physical source graph, with explicit hydrogen atoms.
        path: New CCD destination.
    """
    if path.exists():
        raise FileExistsError(path)
    fresh = Chem.Mol(mol)
    fresh.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261126
    if AllChem.EmbedMolecule(fresh, params) != 0:
        raise ValueError("Freshgraph-onlyCCDconformergenerationfailed")
    lines = [
        f"data_{COMPONENT}",
        "#",
        f"_chem_comp.id {COMPONENT}",
        "_chem_comp.name 'Noncovalent zwitterionic PA6 linear tetramer E-S'",
        "_chem_comp.type non-polymer",
        f"_chem_comp.formula '{rdMolDescriptors.CalcMolFormula(fresh)}'",
        "_chem_comp.mon_nstd_parent_comp_id ?",
        "_chem_comp.pdbx_synonyms ?",
        f"_chem_comp.formula_weight {Descriptors.MolWt(fresh):.4f}",
        "#",
        "loop_",
        "_chem_comp_atom.comp_id",
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.charge",
        "_chem_comp_atom.pdbx_leaving_atom_flag",
        "_chem_comp_atom.pdbx_stereo_config",
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
    ]
    conformer = fresh.GetConformer()
    for atom in fresh.GetAtoms():
        point = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{COMPONENT} {atom.GetProp('_TriposAtomName')} {atom.GetSymbol()} "
            f"{atom.GetFormalCharge()} N N {point.x:.6f} {point.y:.6f} {point.z:.6f}"
        )
    lines += [
        "#",
        "loop_",
        "_chem_comp_bond.comp_id",
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.value_order",
        "_chem_comp_bond.pdbx_aromatic_flag",
    ]
    for bond in fresh.GetBonds():
        order = {1.0: "SING", 2.0: "DOUB"}[bond.GetBondTypeAsDouble()]
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        lines.append(
            f"{COMPONENT} {a.GetProp('_TriposAtomName')} {b.GetProp('_TriposAtomName')} {order} N"
        )
    path.write_text("\n".join([*lines, "#"]) + "\n")


def prepare_graph(output: Path, native_table_path: Path) -> dict:
    """Prepare a graph-only E-S CCD and authoritative per-predictor atom maps.

    Args:
        output: Existing new prediction-preparation directory.
        native_table_path: Actual installed Boltz parser atom-table receipt.

    Returns:
        Graph identity, CCD path and maps for every ligand chain I–L.
    """
    mol = source_molecule(CONFIG / "whole_l4_source.mol2")
    native = json.loads(native_table_path.read_text())
    maps = ligand_maps(mol, native)
    ccd = output / "pa6_l4_noncovalent_es.cif"
    write_es_ccd(mol, ccd)
    result = dict(
        smiles=SMILES,
        component=COMPONENT,
        user_ccd_path=str(ccd.resolve()),
        atoms_total=79,
        heavy_atoms=33,
        formal_charge=0,
        ligand_atom_maps={
            predictor: {chain: mapping for chain in "IJKL"}
            for predictor, mapping in maps.items()
        },
        source_hashes={
            str(p): digest(p)
            for p in [
                CONFIG / "whole_l4_source.mol2",
                CONFIG / "whole_l4_context_mapping.json",
                native_table_path,
                Path(__file__),
                ccd,
            ]
        },
        geometry_input="FreshETKDGgraph-onlyidealconformer;nosourceprotein/TScoordinates",
        covalent_ES_bonds=[],
        amide_selection="Allthreephysicalinternalamidebondsretained;metricsselectseparately",
    )
    (output / "ligand_graph_manifest.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result
