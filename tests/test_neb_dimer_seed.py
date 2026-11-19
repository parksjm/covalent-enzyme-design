"""Focused checks for native NEB validation and local dimer seed geometry."""

import numpy as np
import parmed as pmd
import pytest

from nyl12_covalent.neb_dimer_seed import (
    build_local_dimer_seed,
    final_neb_rows,
    local_heavy_atom_clashes,
    validate_reactive_atom_mapping,
)


def test_local_seed_interpolates_and_normalizes_forward_mode() -> None:
    frames = [
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [8.0, 8.0, 8.0]]),
        np.array([[0.0, 0.2, 0.0], [1.0, 0.4, 0.0], [8.0, 8.0, 8.0]]),
    ]
    seed, mode, norm = build_local_dimer_seed(
        frames, [0, 1], np.array([0, 1]), 1, 2, 0.01
    )
    np.testing.assert_allclose(seed, (frames[0] + frames[1]) / 2)
    assert norm == pytest.approx(np.sqrt(0.20))
    assert np.linalg.norm(mode[[0, 1]] - seed[[0, 1]]) == pytest.approx(0.01)
    assert np.all(mode[0:2, 1] > seed[0:2, 1])
    assert mode[2, 0] == seed[2, 0]


def test_seed_rejects_bad_image_or_rotation_mapping() -> None:
    frames = [np.zeros((2, 3)), np.ones((2, 3))]
    with pytest.raises(ValueError, match="adjacent"):
        build_local_dimer_seed(frames, [0, 1], np.array([0]), 1, 3)
    with pytest.raises(ValueError, match="contained"):
        build_local_dimer_seed(frames, [0], np.array([1]), 1, 2)


def test_final_neb_report_requires_complete_ordered_saves() -> None:
    header = "No climbing image used in cycle 9\nNEB Report\n"
    report = (
        "Img 1 -10.0000000 0.00100 0.00020 1.5\nImg 2 -9.9000000 0.00200 0.00030 1.4\n"
    )
    tail = (
        "Total path length\nSaving img 1 to restart: one\n"
        "Saving img 2 to restart: two\nTesting convergence in cycle 9\n"
    )
    rows = final_neb_rows(header + report + tail, 9, 2)
    assert rows == [(1, -10.0, 0.001, 0.0002), (2, -9.9, 0.002, 0.0003)]
    with pytest.raises(ValueError, match="save loop"):
        final_neb_rows(
            header + report + tail.replace("Saving img 2", "Saved img 2"), 9, 2
        )


def test_reactive_atom_mapping_checks_identity_and_activity() -> None:
    topology = pmd.Structure()
    topology.add_atom(pmd.Atom(name="OG1", atomic_number=8, mass=16), "THR", 929)
    topology.add_atom(pmd.Atom(name="HG1", atomic_number=1, mass=1), "THR", 929)
    declared = {"1": {"name": "OG1", "residue": 1}}
    required = {2: {"name": "HG1", "residue": 1, "atomic_number": 1}}
    assert validate_reactive_atom_mapping(topology, [0, 1], declared, required) == {
        1: 0,
        2: 1,
    }
    with pytest.raises(ValueError, match="identity mismatch"):
        validate_reactive_atom_mapping(
            topology, [0, 1], {"1": {"name": "CA", "residue": 1}}, required
        )


def test_local_clash_report_excludes_reactive_and_bonded_pairs() -> None:
    topology = pmd.Structure()
    for name, number, xyz in [
        ("C", 6, [0.0, 0.0, 0.0]),
        ("N", 7, [1.3, 0.0, 0.0]),
        ("O1", 8, [0.0, 2.0, 0.0]),
        ("O2", 8, [1.4, 2.0, 0.0]),
    ]:
        atom = pmd.Atom(name=name, atomic_number=number, mass=float(number * 2))
        topology.add_atom(atom, "MOL", 1)
    topology.coordinates = np.array(
        [[0.0, 0.0, 0.0], [1.3, 0.0, 0.0], [0.0, 2.0, 0.0], [1.4, 2.0, 0.0]]
    )
    topology.bonds.append(pmd.Bond(topology.atoms[0], topology.atoms[1]))
    contacts = local_heavy_atom_clashes(
        topology, topology.coordinates, {0, 1}, neighborhood_angstrom=4.0
    )
    assert len(contacts["contacts"]) == 1
    assert contacts["contacts"][0]["distance_angstrom"] == pytest.approx(1.4)
    assert contacts["severe_contact_count"] == 1
