"""Synthetic unit tests for TI1 relay-water seed geometry screening."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.relay_water import _contact_limits, _frame, _unit, search_relay_water


def _topology(*, blocker: bool = False, unsupported: bool = False,
              malformed_water: bool = False) -> tuple[SimpleNamespace, np.ndarray]:
    """Build sparse atoms with the residue indices expected by the search."""
    residues = [SimpleNamespace(idx=i, atoms=[]) for i in range(1314)]
    atoms = []
    xyz: list[list[float]] = []

    def add(residue_idx: int, name: str, element: int,
            position: tuple[float, float, float], partners: tuple[int, ...] = ()) -> None:
        residue = residues[residue_idx]
        atom = SimpleNamespace(
            idx=len(atoms), name=name, atomic_number=element, residue=residue, bond_partners=[]
        )
        atoms.append(atom)
        residue.atoms.append(atom)
        xyz.append(list(position))
        atom.bond_partners.extend(atoms[index] for index in partners)

    add(928, "N", 7, (0.0, 0.0, 0.0))
    add(928, "H1", 1, (-0.4, 0.9, 0.0), (0,))
    add(928, "H2", 1, (-0.4, -0.45, 0.78), (0,))
    add(928, "H3", 1, (-0.4, -0.45, -0.78), (0,))
    add(928, "OG1", 8, (2.8, 0.0, 0.0))
    add(1312, "OX", 8, (3.8, 0.0, 0.0))
    if blocker:
        for number, angle in enumerate(np.linspace(0, 2 * np.pi, 12, endpoint=False)):
            add(10, f"C{number}", 6,
                (1.4, 2.48 * np.cos(angle), 2.48 * np.sin(angle)))
    if unsupported:
        add(10, "X", 9, (1.4, 0.0, 3.0))
    water_names = ("O", "H1", "H2") if not malformed_water else ("O", "H1")
    for name in water_names:
        add(1313, name, 8 if name == "O" else 1, (30.0, 30.0, 30.0))
    topology = SimpleNamespace(atoms=atoms, residues=residues)
    return topology, np.asarray(xyz, dtype=float)


def test_search_returns_accepted_water_geometry() -> None:
    """A sparse, unobstructed pocket yields finite accepted water seeds."""
    topology, coordinates = _topology()
    result = search_relay_water(topology, coordinates, "OX", retain=2)
    assert result["local_environment_atom_count"] == 6
    assert result["accepted_oriented_candidate_count"] > 0
    assert len(result["candidates"]) == 2
    candidate = result["candidates"][0]
    assert candidate["accepted"] is True
    oxygen, hydrogen_1, hydrogen_2 = np.asarray(candidate["water_xyz_A"])
    assert np.linalg.norm(oxygen - hydrogen_1) == pytest.approx(0.9572)
    assert np.linalg.norm(oxygen - hydrogen_2) == pytest.approx(0.9572)
    assert candidate["hydrogen_1_donor_target"] in {"OG1", "oxyanion"}
    assert result["water_atom_indices_zero_based"] == [6, 7, 8]


def test_contact_blocker_fails_acceptance_gate() -> None:
    """A steric blocker remains in results but no seed passes the contact gate."""
    topology, coordinates = _topology(blocker=True)
    result = search_relay_water(topology, coordinates, "OX", retain=3)
    assert result["candidates"]
    assert result["accepted_oriented_candidate_count"] == 0
    assert all(candidate["accepted"] is False for candidate in result["candidates"])
    assert result["candidates"][0]["maximum_contact_shortfall_A"] > 0.10


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"malformed_water": True}, "hydrogenated relay water"),
        ({"unsupported": True}, "Unsupported local element"),
    ],
)
def test_invalid_water_inventory_and_unsupported_element_fail(
    kwargs: dict[str, bool], message: str,
) -> None:
    """Invalid explicit water or unsupported local atoms stop the search."""
    topology, coordinates = _topology(**kwargs)
    with pytest.raises(ValueError, match=message):
        search_relay_water(topology, coordinates, "OX")


def test_bad_topology_inventory_and_retain_fail() -> None:
    """The search refuses mismatched inventories and nonpositive candidate counts."""
    topology, coordinates = _topology()
    with pytest.raises(ValueError, match="Unexpected TI1 inventory"):
        search_relay_water(SimpleNamespace(atoms=topology.atoms, residues=[]), coordinates, "OX")
    with pytest.raises(ValueError, match="Unexpected TI1 inventory"):
        search_relay_water(topology, coordinates, "OX", retain=0)


def test_invalid_axis_and_contact_radius_inputs_fail() -> None:
    """Degenerate axes and unsupported elements are rejected directly."""
    with pytest.raises(ValueError, match="coincident atoms"):
        _unit(np.zeros(3))
    with pytest.raises(ValueError, match="coincident atoms"):
        _frame(np.zeros(3))
    topology, _ = _topology(unsupported=True)
    with pytest.raises(ValueError, match="Unsupported local element"):
        _contact_limits(topology, [len(topology.atoms) - 4])
