"""Tests for chemistry-identity PA6 label remapping and PDB normalization."""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Any

import pytest
from rdkit import Chem

from nyl12_covalent import design_component_normalization as normalization

PA6_SMILES = "[NH3+]CCCCCC(=O)NCCCCCC(=O)NCCCCCC(=O)NCCCCCC(=O)[O-]"


def _molecules() -> tuple[Any, Any, dict[str, str], dict[str, str]]:
    molecule = Chem.MolFromSmiles(PA6_SMILES)
    assert molecule is not None
    heavy = Chem.RemoveHs(molecule)
    assert heavy.GetNumAtoms() == 33
    role_indices: dict[str, int] = {}
    for atom in heavy.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygen = [
            bond.GetOtherAtom(atom)
            for bond in atom.GetBonds()
            if bond.GetOtherAtom(atom).GetAtomicNum() == 8
            and bond.GetBondType() == Chem.BondType.DOUBLE
        ]
        nitrogen = [
            bond.GetOtherAtom(atom)
            for bond in atom.GetBonds()
            if bond.GetOtherAtom(atom).GetAtomicNum() == 7
            and bond.GetBondType() == Chem.BondType.SINGLE
        ]
        if oxygen and nitrogen:
            role_indices = {
                "C54": atom.GetIdx(),
                "O49": oxygen[0].GetIdx(),
                "N59": nitrogen[0].GetIdx(),
            }
            break
    assert len(role_indices) == 3
    reference = Chem.Mol(heavy)
    canonical_by_index: dict[int, str] = {}
    other_names = iter(f"X{i:02d}" for i in range(30))
    for atom in reference.GetAtoms():
        role_name = next(
            (role for role, index in role_indices.items() if atom.GetIdx() == index),
            None,
        )
        name = role_name if role_name is not None else next(other_names)
        atom.SetProp("_TriposAtomName", name)
        canonical_by_index[atom.GetIdx()] = name
    order = list(reversed(range(heavy.GetNumAtoms())))
    source = Chem.RenumberAtoms(heavy, order)
    source_to_target: dict[str, str] = {}
    source_elements: dict[str, str] = {}
    for atom in source.GetAtoms():
        source_name = f"B{atom.GetIdx():03d}"
        atom.SetProp("name", source_name)
        old_index = order[atom.GetIdx()]
        source_to_target[source_name] = canonical_by_index[old_index]
        source_elements[source_name] = atom.GetSymbol()
    return source, reference, source_to_target, source_elements


def _pdb_atom(
    record: str,
    serial: int,
    name: str,
    residue: str,
    chain: str,
    number: int,
    element: str,
    xyz: tuple[float, float, float],
) -> str:
    x, y, z = xyz
    return (
        f"{record:<6s}{serial:5d} {name:>4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}"
    )


def _source_pdb(source_elements: dict[str, str]) -> tuple[str, dict[str, str]]:
    lines = []
    proteins: dict[str, str] = {}
    serial = 1
    for chain in "ABCDEFGH":
        line = _pdb_atom("ATOM", serial, "CA", "ALA", chain, 1, "C", (1, 2, 3))
        lines.append(line)
        proteins[chain] = line
        serial += 1
    names = list(source_elements)
    for chain_index, chain in enumerate("IJKL"):
        order = names[chain_index:] + names[:chain_index]
        for index, name in enumerate(order):
            lines.append(
                _pdb_atom(
                    "HETATM",
                    serial,
                    name,
                    "BOL",
                    chain,
                    5,
                    source_elements[name],
                    (float(index), float(chain_index), 0.25),
                )
            )
            serial += 1
    return "\n".join([*lines, "END"]) + "\n", proteins


def test_exact_isomorphism_remaps_labels_and_preserves_amide_roles() -> None:
    source, reference, expected, _ = _molecules()

    result = normalization._match_atom_labels(source, reference)

    assert result["source_to_reference"] == expected
    assert result["graph"] == {
        "heavy_atom_count": 33,
        "heavy_bond_count": 32,
        "heavy_formula": {"C": 24, "N": 4, "O": 5},
        "formal_charge": 0,
        "exact_isomorphism_count": 1,
    }
    assert len(result["bond_mapping"]) == 32
    assert result["role_mapping"]["C54"] in expected
    assert result["role_mapping"]["O49"] in expected
    assert result["role_mapping"]["N59"] in expected


def test_formal_charge_mismatch_is_rejected() -> None:
    source, reference, _, _ = _molecules()
    source.GetAtomWithIdx(0).SetFormalCharge(1)
    with pytest.raises(ValueError, match="formal charges must be zero"):
        normalization._match_atom_labels(source, reference)


def test_ambiguous_exact_graph_mappings_are_rejected() -> None:
    source, reference, _, _ = _molecules()
    candidates = reference.GetSubstructMatches(source, uniquify=False)
    assert len(candidates) == 1
    with pytest.raises(ValueError, match="found 2"):
        normalization._unique_graph_match(
            source, reference, (candidates[0], candidates[0])
        )


def test_duplicate_source_atom_names_are_rejected() -> None:
    source, reference, _, _ = _molecules()
    source.GetAtomWithIdx(1).SetProp("name", source.GetAtomWithIdx(0).GetProp("name"))
    with pytest.raises(ValueError, match="non-empty and unique"):
        normalization._match_atom_labels(source, reference)


def test_pdb_normalization_uses_graph_map_and_preserves_protein_and_xyz(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source_molecule, reference_molecule, expected, source_elements = _molecules()
    text, protein_lines = _source_pdb(source_elements)
    pdb = tmp_path / "boltz.pdb"
    pickle_path = tmp_path / "component.pkl"
    mol2 = tmp_path / "reference.mol2"
    pdb.write_text(text)
    pickle_path.write_bytes(pickle.dumps({"trusted fixture": True}))
    mol2.write_text("reference fixture\n")
    monkeypatch.setattr(
        normalization, "_load_component_molecule", lambda _: source_molecule
    )
    monkeypatch.setattr(
        normalization, "_load_reference_molecule", lambda _: reference_molecule
    )
    output_dir = tmp_path / "normalized"

    result = normalization.normalize_model(pdb, pickle_path, mol2, output_dir)

    normalized = (output_dir / "normalized.pdb").read_text().splitlines()
    assert {line for line in normalized if line.startswith("ATOM  ")} == set(
        protein_lines.values()
    )
    ligands = [line for line in normalized if line.startswith("HETATM")]
    assert len(ligands) == 132
    original_by_serial = {
        (line[21], int(line[6:11])): line
        for line in text.splitlines()
        if line.startswith("HETATM")
    }
    map_by_serial = {
        (row["source_chain"], row["source_serial"]): row["target_atom"]
        for row in result["pdb_atom_mapping"]
    }
    for line in ligands:
        assert line[17:20] == "LIG"
        assert line[22:26] == "   1"
        key = (line[21], int(line[6:11]))
        assert line[12:16].strip() == map_by_serial[key]
        assert line[30:54] == original_by_serial[key][30:54]
        assert line[12:16].strip() in expected.values()
    assert result["role_mapping"]["amide_bonds"] == {
        "C54-O49": "DOUBLE",
        "C54-N59": "SINGLE",
    }
    assert len(result["pdb_atom_mapping"]) == 132
    assert (
        json.loads((output_dir / "mapping.json").read_text())["graph"][
            "heavy_atom_count"
        ]
        == 33
    )


@pytest.mark.parametrize(
    "bad_field, bad_value, message",
    [
        ("altloc", "A", "Alternate locations"),
        ("coordinate", "     nan", "non-finite"),
    ],
)
def test_bad_pdb_features_are_rejected(
    bad_field: str, bad_value: str, message: str
) -> None:
    _, _, _, elements = _molecules()
    lines = _source_pdb(elements)[0].splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("HETATM"))
    line = lines[index]
    lines[index] = (
        f"{line[:16]}{bad_value}{line[17:]}"
        if bad_field == "altloc"
        else f"{line[:30]}{bad_value}{line[38:]}"
    )
    with pytest.raises(ValueError, match=message):
        normalization._normalize_pdb_lines(
            lines, {name: name for name in elements}, elements
        )


def test_pdb_component_must_match_pickle_atom_label_set() -> None:
    _, _, _, elements = _molecules()
    lines = _source_pdb(elements)[0].splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("HETATM"))
    lines[index] = f"{lines[index][:12]}MISS{lines[index][16:]}"
    with pytest.raises(ValueError, match="source atom labels differ"):
        normalization._normalize_pdb_lines(
            lines, {name: name for name in elements}, elements
        )
