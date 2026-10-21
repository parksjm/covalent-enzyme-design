"""Focused tests for pure rigid geometry and relay repair screens."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.relay_repair_geometry import (
    axis_local_environment_indices,
    component_after_cut,
    donor_geometry,
    graph_exclusions,
    moved_atom_contact_screen,
    oxygen_grid,
    pair_contact_limit,
    rotate_points,
    unit,
    water_geometry,
)


def _graph() -> SimpleNamespace:
    """Create a small bonded graph with representative local elements."""
    atoms: list[SimpleNamespace] = []

    def add(element: int, parent: int | None = None) -> int:
        atom = SimpleNamespace(idx=len(atoms), atomic_number=element, bond_partners=[])
        atoms.append(atom)
        if parent is not None:
            atom.bond_partners.append(atoms[parent])
            atoms[parent].bond_partners.append(atom)
        return atom.idx

    ca = add(6)
    cb = add(6, ca)
    cg = add(6, cb)
    h_nonpolar = add(1, cg)
    n = add(7)
    h_polar = add(1, n)
    o = add(8)
    sulfur = add(16)
    isolated_c = add(6)
    return SimpleNamespace(
        atoms=atoms,
        ca=ca,
        cb=cb,
        cg=cg,
        h_nonpolar=h_nonpolar,
        n=n,
        h_polar=h_polar,
        o=o,
        sulfur=sulfur,
        isolated_c=isolated_c,
    )


def test_unit_rotation_and_graph_cut_preserve_rigid_distances() -> None:
    """A cut component rotates as one unit without changing internal lengths."""
    graph = _graph()
    coords = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0], [4.0, 0.0, 0.0], [4.0, 1.0, 0.0],
        [5.0, 0.0, 0.0], [5.0, 1.0, 0.0], [8.0, 0.0, 0.0],
    ])
    selected = component_after_cut(graph.atoms[graph.cb], graph.atoms[graph.ca])
    assert selected == [graph.cb, graph.cg, graph.h_nonpolar]
    assert graph_exclusions(graph.atoms[graph.ca]) == {graph.ca, graph.cb, graph.cg}
    rotated = rotate_points(coords, selected, coords[graph.ca], [1, 0, 0], 90)
    assert np.array_equal(rotated[graph.ca], coords[graph.ca])
    assert np.array_equal(rotated[graph.isolated_c], coords[graph.isolated_c])
    assert np.linalg.norm(rotated[graph.cb] - rotated[graph.cg]) == pytest.approx(
        np.linalg.norm(coords[graph.cb] - coords[graph.cg])
    )
    assert unit(np.asarray([0.0, 3.0, 4.0])) == pytest.approx([0.0, 0.6, 0.8])
    with pytest.raises(ValueError, match="coincident"):
        unit(np.zeros(3))


@pytest.mark.parametrize(
    ("left", "expected"),
    [
        ("H-H", 1.40),
        ("H-C", 2.10),
        ("H-polar", 1.55),
        ("H-nonpolar", 2.00),
        ("H-S", 2.30),
        ("C-C", 2.80),
        ("C-polar", 2.75),
        ("C-S", 3.00),
        ("polar-polar", 2.50),
        ("polar-S", 2.95),
        ("S-S", 3.20),
    ],
)
def test_pair_contact_limits_are_explicit_and_symmetric(
    left: str, expected: float
) -> None:
    """Each requested atom-pair screening threshold is applied symmetrically."""
    graph = _graph()
    names = {
        "H-H": (graph.atoms[graph.h_polar], graph.atoms[graph.h_nonpolar]),
        "H-C": (graph.atoms[graph.h_polar], graph.atoms[graph.ca]),
        "H-polar": (graph.atoms[graph.h_polar], graph.atoms[graph.o]),
        "H-nonpolar": (graph.atoms[graph.h_nonpolar], graph.atoms[graph.n]),
        "H-S": (graph.atoms[graph.h_polar], graph.atoms[graph.sulfur]),
        "C-C": (graph.atoms[graph.ca], graph.atoms[graph.isolated_c]),
        "C-polar": (graph.atoms[graph.ca], graph.atoms[graph.o]),
        "C-S": (graph.atoms[graph.ca], graph.atoms[graph.sulfur]),
        "polar-polar": (graph.atoms[graph.n], graph.atoms[graph.o]),
        "polar-S": (graph.atoms[graph.o], graph.atoms[graph.sulfur]),
        "S-S": (graph.atoms[graph.sulfur], graph.atoms[graph.sulfur]),
    }
    a, b = names[left]
    assert pair_contact_limit(a, b) == pytest.approx(expected)
    assert pair_contact_limit(b, a) == pytest.approx(expected)


def test_unsupported_element_and_graph_exclusions_are_fail_closed() -> None:
    """Unknown elements raise; contact screen excludes only graph distance 1-2/1-3."""
    graph = _graph()
    bad = SimpleNamespace(atomic_number=9, bond_partners=[])
    with pytest.raises(ValueError, match="Unsupported"):
        pair_contact_limit(graph.atoms[graph.ca], bad)
    graph.atoms.append(SimpleNamespace(
        idx=len(graph.atoms), atomic_number=6, bond_partners=[]
    ))
    env = graph.isolated_c
    coords = np.zeros((len(graph.atoms), 3), dtype=float)
    coords[graph.isolated_c] = [2.5, 0.0, 0.0]
    result = moved_atom_contact_screen(
        graph, coords, [graph.ca], [graph.ca, graph.cb, graph.cg, env]
    )
    assert result["accepted"] is False
    assert result["violation_count"] == 1
    assert result["violations"][0]["environment_atom_index_zero_based"] == env


def test_axis_local_environment_bound_keeps_possible_contacts_only() -> None:
    """The filter includes the exact radial motion plus maximum-threshold bound."""
    coords = np.asarray([[0.0, 0.0, 0.0], [0.0, 7.0, 0.0], [0.0, 7.21, 0.0]])
    selected = axis_local_environment_indices(
        coords, [0, 1, 2], np.zeros(3), np.asarray([1.0, 0.0, 0.0]), 4.0
    )
    assert selected == [0, 1]


def test_oxygen_grid_water_geometry_and_donor_angle() -> None:
    """Grid radii, rigid water geometry, and angle convention are explicit."""
    nitrogen = np.asarray([0.0, 0.0, 0.0])
    og1 = np.asarray([2.8, 0.0, 0.0])
    sites = oxygen_grid(nitrogen, og1)
    assert sites
    oxygen, n_distance, og_distance, *_ = sites[0]
    assert np.linalg.norm(oxygen - nitrogen) == pytest.approx(n_distance)
    assert np.linalg.norm(oxygen - og1) == pytest.approx(og_distance)
    water = water_geometry(oxygen, og1, 120)
    assert np.linalg.norm(water[1] - water[0]) == pytest.approx(0.9572)
    assert np.linalg.norm(water[2] - water[0]) == pytest.approx(0.9572)
    oh1 = water[1] - water[0]
    oh2 = water[2] - water[0]
    theta = np.degrees(np.arccos(np.dot(oh1, oh2) / np.linalg.norm(oh1) / np.linalg.norm(oh2)))
    assert theta == pytest.approx(104.52)
    distance, angle = donor_geometry([0, 0, 0], [1, 0, 0], [2.8, 0, 0])
    assert distance == pytest.approx(1.8)
    assert angle == pytest.approx(180.0)
    _, wrong_angle = donor_geometry([0, 0, 0], [1, 0, 0], [0.2, 0, 0])
    assert wrong_angle == pytest.approx(0.0)
