"""Synthetic Amber topology coverage for the read-only ES auditor."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import parmed
import pytest
from parmed.topologyobjects import Atom, AtomType, Bond, BondType

from nyl12_covalent.es_topology_audit import audit_es_topology


@pytest.fixture
def amber_fixture(tmp_path: Path) -> tuple[Path, dict[str, Any]]:
    """Write a minimal actual Amber file with explicit hydrogen bonds."""
    structure = parmed.Structure()
    specification = [
        ("N", 7, "THR", 1),
        ("H1", 1, "THR", 1),
        ("H2", 1, "THR", 1),
        ("OG1", 8, "THR", 1),
        ("HG1", 1, "THR", 1),
        ("C", 6, "LIG", 2),
        ("O", 8, "LIG", 2),
        ("N", 7, "LIG", 2),
        ("H", 1, "LIG", 2),
        ("N", 7, "ASN", 3),
        ("H", 1, "ASN", 3),
        ("ND2", 7, "ASN", 3),
        ("HD21", 1, "ASN", 3),
    ]
    masses = {1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999}
    for index, (name, element, residue, number) in enumerate(specification):
        atom_type = AtomType(f"T{element}", element, masses[element], element)
        atom_type.set_lj_params(0.1, 1.5)
        atom = Atom(
            name=name,
            type=f"T{element}",
            atomic_number=element,
            mass=masses[element],
            charge=0.125 if index == 0 else 0,
        )
        atom.atom_type = atom_type
        structure.add_atom(atom, residue, number)
    bond_type = BondType(100.0, 1.0)
    structure.bond_types.append(bond_type)
    for left, right in [
        (0, 1),
        (0, 2),
        (0, 3),
        (3, 4),
        (5, 6),
        (5, 7),
        (7, 8),
        (9, 10),
        (9, 11),
        (11, 12),
    ]:
        structure.bonds.append(
            Bond(structure.atoms[left], structure.atoms[right], type=bond_type)
        )
    path = tmp_path / "es.prmtop"
    parmed.amber.AmberParm.from_structure(structure).save(str(path))
    mapping = {
        "nucleophile": {"residue": 1, "resname": "THR", "atom": "OG1"},
        "scissile_carbon": {"residue": 2, "resname": "LIG", "atom": "C"},
        "scissile_oxygen": {"residue": 2, "resname": "LIG", "atom": "O"},
        "leaving_nitrogen": {"residue": 2, "resname": "LIG", "atom": "N"},
        "oxyanion_donor_backbone": {"residue": 3, "resname": "ASN", "atom": "N"},
        "oxyanion_donor_sidechain": {"residue": 3, "resname": "ASN", "atom": "ND2"},
        "thr_n_terminus": {"residue": 1, "resname": "THR", "atom": "N"},
    }
    return path, mapping


def test_actual_amber_inventory(amber_fixture: tuple[Path, dict[str, Any]]) -> None:
    path, mapping = amber_fixture
    before = path.read_bytes()
    result = audit_es_topology(path, mapping, 2)
    assert path.read_bytes() == before
    assert result["selected_atoms"]["nucleophile"]["index"] == 4
    assert result["selected_atoms"]["nucleophile"]["explicit_hydrogen_neighbors"] == [5]
    assert result["protonation_observations"]["thr_n_terminus"] == {
        "explicit_hydrogen_count": 2
    }
    assert result["mapped_residues"][0]["net_charge"] == pytest.approx(0.125)
    assert result["ligand"]["connected_components"] == [[6, 7, 8, 9]]
    assert result["ligand"]["bonds"] == [[6, 7], [6, 8], [8, 9]]


@pytest.mark.parametrize(
    "mutation",
    ["missing", "extra", "bad_row", "residue", "bool", "resname", "atom", "empty_atom"],
)
def test_mapping_errors(
    amber_fixture: tuple[Path, dict[str, Any]],
    mutation: str,
) -> None:
    path, original = amber_fixture
    mapping = copy.deepcopy(original)
    if mutation == "missing":
        del mapping["nucleophile"]
    elif mutation == "extra":
        mapping["unknown"] = mapping["nucleophile"]
    elif mutation == "bad_row":
        mapping["nucleophile"]["unknown"] = 1
    else:
        key, value = {
            "residue": ("residue", 99),
            "bool": ("residue", True),
            "resname": ("resname", "SER"),
            "atom": ("atom", "NOPE"),
            "empty_atom": ("atom", ""),
        }[mutation]
        mapping["nucleophile"][key] = value
    with pytest.raises(ValueError):
        audit_es_topology(path, mapping, 2)


@pytest.mark.parametrize(
    "mutation",
    [
        "covalent",
        "disconnected",
        "missing_bond",
        "external",
        "duplicate",
        "nonfinite",
        "missing_h",
    ],
)
def test_corrupt_topology(
    amber_fixture: tuple[Path, dict[str, Any]],
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    path, mapping = amber_fixture
    structure = parmed.amber.AmberParm(str(path))
    if mutation in {"covalent", "external"}:
        left = 3 if mutation == "covalent" else 9
        structure.bonds.append(Bond(structure.atoms[left], structure.atoms[5]))
    elif mutation in {"disconnected", "missing_bond", "missing_h"}:
        indices = {"disconnected": (7, 8), "missing_bond": (5, 6), "missing_h": (3, 4)}[
            mutation
        ]
        bond = next(
            bond
            for bond in structure.bonds
            if {bond.atom1.idx, bond.atom2.idx} == set(indices)
        )
        bond.delete()
        structure.bonds.remove(bond)
    elif mutation == "duplicate":
        structure.atoms[4].name = "OG1"
    else:
        structure.atoms[3].charge = float("nan")
    monkeypatch.setattr(parmed.amber, "AmberParm", lambda _: structure)
    if mutation == "missing_h":
        result = audit_es_topology(path, mapping, 2)
        assert result["protonation_observations"]["nucleophile"] == {
            "explicit_hydrogen_count": 0
        }
    else:
        with pytest.raises(ValueError):
            audit_es_topology(path, mapping, 2)


def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        audit_es_topology(tmp_path / "absent.prmtop", {}, 2)


@pytest.mark.parametrize(
    "label",
    [
        "nucleophile",
        "scissile_carbon",
        "scissile_oxygen",
        "leaving_nitrogen",
        "oxyanion_donor_backbone",
        "oxyanion_donor_sidechain",
        "thr_n_terminus",
    ],
)
def test_wrong_role_element(
    amber_fixture: tuple[Path, dict[str, Any]],
    monkeypatch: pytest.MonkeyPatch,
    label: str,
) -> None:
    path, mapping = amber_fixture
    structure = parmed.amber.AmberParm(str(path))
    row = mapping[label]
    atom = next(
        atom
        for atom in structure.residues[row["residue"] - 1].atoms
        if atom.name == row["atom"]
    )
    # The leaving group remains adjacent to scissile C, but is now carbon.
    atom.atomic_number = 6 if atom.atomic_number != 6 else 7
    monkeypatch.setattr(parmed.amber, "AmberParm", lambda _: structure)
    with pytest.raises(ValueError, match=f"Unexpected element for {label}"):
        audit_es_topology(path, mapping, 2)


def test_duplicate_roles(amber_fixture: tuple[Path, dict[str, Any]]) -> None:
    path, mapping = amber_fixture
    mapping["oxyanion_donor_sidechain"] = mapping["oxyanion_donor_backbone"].copy()
    with pytest.raises(ValueError, match="Distinct roles"):
        audit_es_topology(path, mapping, 2)


def test_adjacent_carbon_cannot_be_leaving_nitrogen(
    amber_fixture: tuple[Path, dict[str, Any]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path, mapping = amber_fixture
    structure = parmed.amber.AmberParm(str(path))
    adjacent = structure.atoms[8]
    adjacent.atomic_number = 6
    adjacent.name = "CA"
    structure.bonds.append(Bond(structure.atoms[5], adjacent))
    mapping["leaving_nitrogen"]["atom"] = "CA"
    monkeypatch.setattr(parmed.amber, "AmberParm", lambda _: structure)
    with pytest.raises(ValueError, match="Unexpected element for leaving_nitrogen"):
        audit_es_topology(path, mapping, 2)


def test_terminal_nitrogen_wrong_residue(
    amber_fixture: tuple[Path, dict[str, Any]],
) -> None:
    path, mapping = amber_fixture
    mapping["thr_n_terminus"] = {"residue": 3, "resname": "ASN", "atom": "ND2"}
    mapping["oxyanion_donor_sidechain"] = {
        "residue": 1,
        "resname": "THR",
        "atom": "N",
    }
    with pytest.raises(ValueError, match="must share a residue"):
        audit_es_topology(path, mapping, 2)
