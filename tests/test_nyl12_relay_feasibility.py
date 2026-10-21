"""Synthetic tests for the bounded relay oxygen feasibility diagnostic."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.relay_feasibility import (
    _validate_inventory,
    rank_oxygen_sites,
)


def _topology() -> tuple[SimpleNamespace, np.ndarray]:
    """Build a small graph embedded in the validated 1314-residue inventory."""
    residues = [SimpleNamespace(idx=i, name="ALA", atoms=[]) for i in range(1314)]
    residues[185].name = "ASN"
    residues[928].name = "THR"
    residues[1313].name = "WAT"
    atoms: list[SimpleNamespace] = []
    coordinates: list[tuple[float, float, float]] = []

    def add(
        residue_index: int,
        name: str,
        element: int,
        position: tuple[float, float, float],
        parent: int | None = None,
    ) -> SimpleNamespace:
        residue = residues[residue_index]
        atom = SimpleNamespace(
            idx=len(atoms),
            name=name,
            atomic_number=element,
            residue=residue,
            bond_partners=[],
        )
        atoms.append(atom)
        residue.atoms.append(atom)
        coordinates.append(position)
        if parent is not None:
            atom.bond_partners.append(atoms[parent])
            atoms[parent].bond_partners.append(atom)
        return atom

    # ASN186 backbone and sidechain, including sidechain hydrogens.
    add(185, "N", 7, (50.0, 0.0, 0.0))
    ca = add(185, "CA", 6, (51.0, 0.0, 0.0))
    add(185, "C", 6, (52.0, 0.0, 0.0), ca.idx)
    add(185, "H", 1, (50.0, 1.0, 0.0))
    cb = add(185, "CB", 6, (51.0, 1.0, 0.0), ca.idx)
    cg = add(185, "CG", 6, (51.0, 2.0, 0.0), cb.idx)
    od1 = add(185, "OD1", 8, (52.0, 2.0, 0.0), cg.idx)
    nd2 = add(185, "ND2", 7, (50.0, 2.0, 0.0), cg.idx)
    add(185, "HD21", 1, (50.0, 3.0, 0.0), nd2.idx)
    add(185, "HD22", 1, (49.0, 2.0, 0.0), nd2.idx)
    add(185, "HD", 1, (51.0, 3.0, 0.0), od1.idx)

    # Reactive Thr929 with NH3+ hydrogens.
    nitrogen = add(928, "N", 7, (0.0, 0.0, 0.0))
    add(928, "CA", 6, (-1.0, 0.0, 0.0), nitrogen.idx)
    add(928, "OG1", 8, (2.8, 0.0, 0.0))
    for name, position in (
        ("H1", (-0.4, 0.9, 0.0)),
        ("H2", (-0.4, -0.45, 0.78)),
        ("H3", (-0.4, -0.45, -0.78)),
    ):
        add(928, name, 1, position, nitrogen.idx)
    add(1313, "O", 8, (30.0, 30.0, 30.0))
    add(1313, "H1", 1, (30.9, 30.0, 30.0))
    add(1313, "H2", 1, (29.7, 30.9, 30.0))
    return SimpleNamespace(atoms=atoms, residues=residues), np.asarray(coordinates)


def test_fixed_grid_reports_eligible_sites_and_spatially_distinct_ranking() -> None:
    """Clear environments admit sites and ranking has no duplicate points."""
    result = rank_oxygen_sites(
        np.array([0.0, 0.0, 0.0]),
        np.array([2.8, 0.0, 0.0]),
        np.empty((0, 3)),
        np.empty(0),
        retain=10,
    )
    assert result["grid_site_count"] > 0
    assert result["eligible_site_count"] == result["grid_site_count"]
    assert len(result["sites"]) == 10
    points = [np.asarray(site["oxygen_xyz_A"]) for site in result["sites"]]
    assert all(
        np.linalg.norm(a - b) >= 0.4
        for i, a in enumerate(points)
        for b in points[i + 1 :]
    )


def test_fixed_grid_shortfall_is_positive_and_ranked() -> None:
    """A close carbon blocks sites and reports its positive distance deficit."""
    result = rank_oxygen_sites(
        np.array([0.0, 0.0, 0.0]),
        np.array([2.8, 0.0, 0.0]),
        np.array([[1.4, 0.0, 2.4]]),
        np.array([100.0]),
    )
    shortfalls = [site["maximum_positive_shortfall_A"] for site in result["sites"]]
    assert shortfalls == sorted(shortfalls)
    assert result["best_max_shortfall_A"] > 0.10
    assert result["eligible_site_count"] < result["grid_site_count"]


def test_optimistic_removal_selection_follows_bonds_from_asn_cb_only() -> None:
    """ASN backbone atoms remain while sidechain and beta-N hydrogens are selected."""
    topology, coordinates = _topology()
    inventory, removed, waters = _validate_inventory(topology, coordinates)
    labels = {(topology.atoms[i].residue.idx, topology.atoms[i].name) for i in removed}
    assert inventory["finite_coordinates"] is True
    assert (185, "CB") in labels and (185, "HD21") in labels and (185, "HD22") in labels
    assert (
        (185, "CA") not in labels
        and (185, "N") not in labels
        and (185, "H") not in labels
    )
    assert {(928, "H1"), (928, "H2"), (928, "H3")} <= labels
    assert len(waters) == 3


@pytest.mark.parametrize("bad_coordinate", [np.nan, np.inf])
def test_inventory_rejects_nonfinite_coordinates(bad_coordinate: float) -> None:
    """Invalid coordinate arrays fail before any steric calculation."""
    topology, coordinates = _topology()
    coordinates[0, 0] = bad_coordinate
    with pytest.raises(ValueError, match="finite"):
        _validate_inventory(topology, coordinates)
