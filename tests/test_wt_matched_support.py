"""Meaningful coordinate/PBC/graph checks for matched state preparation."""

from pathlib import Path

import numpy as np
import parmed as pmd
import pytest

from nyl12_covalent.wt_matched_support import (
    bond_keys,
    digest,
    distant_sodium,
    write_full_pdb,
)


def system() -> pmd.Structure:
    """Make a small periodic solute plus five ions and a water."""
    s = pmd.Structure()
    for name, residue, x in [
        ("C", "ALA", 0),
        ("Na+", "Na+", 9.9),
        ("Na+", "Na+", 2),
        ("Na+", "Na+", 3),
        ("Na+", "Na+", 4),
        ("Na+", "Na+", 5),
        ("O", "WAT", 7),
    ]:
        a = pmd.Atom(name=name, atomic_number=6 if name == "C" else 11)
        s.add_atom(a, residue, len(s.residues) + 1)
        a.xx, a.xy, a.xz = x, 0, 0
    s.box = [10, 10, 10, 90, 90, 90]
    return s


def test_periodic_ion_ranking_uses_minimum_image() -> None:
    s = system()
    # Ion9.9 is only0.1A away across boundary, so retain it.
    assert distant_sodium(s) == [6, 5, 4, 3]


def test_full_pdb_retains_water_and_mapping(tmp_path: Path) -> None:
    s = system()
    p = tmp_path / "full.pdb"
    mapping = write_full_pdb(s, p, [3], False)
    assert mapping == [1, 2, 4, 5, 6, 7]
    assert "WAT" in p.read_text()
    assert len([x for x in p.read_text().splitlines() if x.startswith("ATOM")]) == 6
    assert len(digest(p)) == 64


def test_not_enough_ions_fails() -> None:
    s = system()
    s.strip(":3-6")
    with pytest.raises(ValueError, match="four sodium"):
        distant_sodium(s)


def test_bond_keys_use_parent_residues() -> None:
    s = system()
    s.bonds.append(pmd.Bond(s.atoms[0], s.atoms[1]))
    assert bond_keys(s, [101, 102, 103, 104, 105, 106, 107]) == {
        ((101, "C"), (102, "Na+"))
    }


def test_restore_shared_coordinates_and_graph(monkeypatch: pytest.MonkeyPatch) -> None:
    from nyl12_covalent import wt_matched_support as support

    monkeypatch.setattr(support, "BETA", ())
    s = pmd.Structure()
    for name, charge, x in [("N", -8.0, 0.0), ("H1", 0.0, 1.0)]:
        atom = pmd.Atom(name=name, atomic_number=7 if name == "N" else 1, charge=charge)
        s.add_atom(atom, "ALA", 1)
        atom.xx, atom.xy, atom.xz = x, 0, 0
    s.bonds.append(pmd.Bond(*s.atoms))
    s.box = [10, 10, 10, 90, 90, 90]
    product = s.copy(pmd.Structure)
    product.coordinates = s.coordinates + 2
    report = support.restore_and_check(s, product, [1], False, False)
    assert np.array_equal(product.coordinates, s.coordinates)
    assert report["dry_charge"] == -8
    assert report["shared_bond_graph_exact"]
    product.bonds.clear()
    with pytest.raises(ValueError, match="bond graph"):
        support.restore_and_check(s, product, [1], False, False)


def test_new_proton_graph_and_superposition(monkeypatch: pytest.MonkeyPatch) -> None:
    from nyl12_covalent import wt_matched_support as support

    monkeypatch.setattr(support, "BETA", ())
    monkeypatch.setattr(support, "ASH", (1,))
    parent = pmd.Structure()
    atom = pmd.Atom(name="OD2", atomic_number=8, charge=-4)
    parent.add_atom(atom, "ASP", 1)
    parent.coordinates = [[0, 0, 0]]
    parent.box = [10, 10, 10, 90, 90, 90]
    product = parent.copy(pmd.Structure)
    product.residues[0].name = "ASH"
    h = pmd.Atom(name="HD2", atomic_number=1, charge=0)
    product.add_atom(h, "ASH", 1)
    product.coordinates = [[3, 0, 0], [1, 0, 0]]
    product.bonds.append(pmd.Bond(product.atoms[0], h))
    report = support.restore_and_check(parent, product, [1], False, True)
    assert report["new_hydrogens"] == [(1, "HD2")]
    assert report["new_proton_min_distance_A"] == pytest.approx(1)
    product.coordinates = [[0, 0, 0], [0.1, 0, 0]]
    with pytest.raises(ValueError, match="superposition"):
        support.restore_and_check(parent, product, [1], False, True)
