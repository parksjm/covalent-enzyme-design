"""Focused tests for identity-safe ASH272 QM-region charge audits."""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from nyl12_covalent.ash272_qmmm import (
    audit_qm_region,
    compact_mask,
    compact_residue_mask,
    remap_active_indices,
)


@dataclass
class Residue:
    """Minimal topology residue for QM-region audit tests."""

    idx: int
    name: str
    atoms: list[Atom] = field(default_factory=list)


@dataclass
class Atom:
    """Minimal topology atom for QM-region audit tests."""

    idx: int
    name: str
    element_name: str
    residue: Residue


@dataclass
class Bond:
    """Minimal topology bond for QM-region audit tests."""

    atom1: Atom
    atom2: Atom


@dataclass
class Topology:
    """Minimal topology fixture."""

    residues: list[Residue]
    atoms: list[Atom]
    bonds: list[Bond]


def _topology(protonated: bool) -> Topology:
    """Build a two-residue Asp/Ash topology with one QM/MM link cut."""
    cap = Residue(0, "ALA")
    site = Residue(1, "ASH" if protonated else "ASP")
    atoms: list[Atom] = []

    def add(residue: Residue, name: str, element: str) -> Atom:
        atom = Atom(len(atoms), name, element, residue)
        residue.atoms.append(atom)
        atoms.append(atom)
        return atom

    add(cap, "CA", "C")
    ca = add(site, "CA", "C")
    cb = add(site, "CB", "C")
    hb2 = add(site, "HB2", "H")
    hb3 = add(site, "HB3", "H")
    cg = add(site, "CG", "C")
    od1 = add(site, "OD1", "O")
    od2 = add(site, "OD2", "O")
    hd2 = add(site, "HD2", "H") if protonated else None
    bonds = [
        Bond(ca, cb),
        Bond(cb, hb2),
        Bond(cb, hb3),
        Bond(cb, cg),
        Bond(cg, od1),
        Bond(cg, od2),
    ]
    if hd2 is not None:
        bonds.append(Bond(od2, hd2))
    return Topology([cap, site], atoms, bonds)


def test_ash_charge_includes_hd2_and_preserves_electron_count() -> None:
    """An included Asp-to-Ash proton adds one nucleus and one formal charge."""
    old, new = _topology(False), _topology(True)
    old_qm = [atom.idx for atom in old.residues[1].atoms if atom.name != "CA"]
    baseline = audit_qm_region(
        old, old, old_qm, {}, reference_charge=-1, expected_link_h_count=1
    )
    audit = audit_qm_region(
        old,
        new,
        old_qm,
        {2: "ASH"},
        reference_charge=-1,
        expected_link_h_count=1,
    )

    assert audit.real_atom_count == len(old_qm) + 1
    assert audit.nuclear_charge_with_links == baseline.nuclear_charge_with_links + 1
    assert audit.electron_count == baseline.electron_count
    assert audit.formal_charge == 0
    assert audit.link_h_count == 1
    assert audit.electron_count == audit.nuclear_charge_with_links
    assert audit.multiplicity == 1
    assert new.residues[1].atoms[-1].idx in audit.atom_indices


def test_protonation_outside_qm_does_not_change_qm_charge() -> None:
    """A declared Ash in the MM environment contributes no QM formal delta."""
    old, new = _topology(False), _topology(True)
    old_qm = [atom.idx for atom in old.residues[0].atoms]
    audit = audit_qm_region(
        old,
        new,
        old_qm,
        {2: "ASH"},
        reference_charge=-2,
        expected_link_h_count=0,
    )

    assert audit.formal_charge == -2
    assert audit.real_atom_count == 1
    assert audit.link_h_count == 0


def test_active_mask_adds_hd2_when_od2_is_mobile() -> None:
    """The new ASH proton follows its selected parent oxygen."""
    old, new = _topology(False), _topology(True)
    od2_index = next(atom.idx for atom in old.residues[1].atoms if atom.name == "OD2")
    remapped = remap_active_indices(old, new, [od2_index], {2: "ASH"})

    hd2_index = next(atom.idx for atom in new.residues[1].atoms if atom.name == "HD2")
    od2_new_index = next(
        atom.idx for atom in new.residues[1].atoms if atom.name == "OD2"
    )
    assert set(remapped) == {hd2_index, od2_new_index}


def test_declared_protonation_must_match_topology() -> None:
    """A mislabeled QM state fails before preparing an input file."""
    old, new = _topology(False), _topology(False)
    old_qm = [atom.idx for atom in old.residues[1].atoms if atom.name != "CA"]

    with pytest.raises(ValueError, match="expected declared ASH"):
        audit_qm_region(
            old,
            new,
            old_qm,
            {2: "ASH"},
            reference_charge=-1,
            expected_link_h_count=1,
        )


def test_compact_mask_merges_contiguous_indices() -> None:
    """The generated mask is deterministic and one-based."""
    assert compact_mask([5, 1, 2, 3, 9, 9]) == "@1-3,5,9"


def test_compact_residue_mask_preserves_complete_residues() -> None:
    """Complete residue selections can use shorter exact Amber masks."""
    topology = _topology(False)
    selected = [atom.idx for atom in topology.residues[1].atoms]

    assert compact_residue_mask(topology, selected) == ":2"


def test_compact_residue_mask_rejects_partial_residues() -> None:
    """A residue expression must not widen the movable atom selection."""
    topology = _topology(False)

    with pytest.raises(ValueError, match="partial residue"):
        compact_residue_mask(topology, [topology.residues[1].atoms[0].idx])


def test_mask_rejects_empty_or_nonpositive_indices() -> None:
    """Invalid masks fail explicitly instead of creating a fallback selection."""
    with pytest.raises(ValueError, match="positive one-based"):
        compact_mask([0, 1])
    with pytest.raises(ValueError, match="requires positive"):
        compact_mask([])
