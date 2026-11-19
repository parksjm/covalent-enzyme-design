"""Verify local acylation guesses preserve chemical and movement inventories."""

from dataclasses import replace

import numpy as np
import parmed
import pytest

from test_acylenzyme_seed import seed  # noqa: F401
from nyl12_covalent.relay_acyl_seed import (
    RelayAcylAtomMap,
    build_balanced_acyl_seed,
)


@pytest.fixture
def relay(seed: tuple) -> tuple:  # noqa: F811
    topology, xyz, acyl = seed
    xyz[4] = xyz[3] + [0, 0, 1]
    atoms = [
        ("O", 8),
        ("H1", 1),
        ("H2", 1),
        ("OD2", 8),
        ("HD2", 1),
        ("OD1", 8),
        ("OD2", 8),
    ]
    for i, (name, number) in enumerate(atoms):
        topology.add_atom(parmed.Atom(name=name, atomic_number=number), "WAT", 20 + i)
    for a, b in ((14, 15), (14, 16), (17, 18)):
        topology.bonds.append(parmed.Bond(topology.atoms[a], topology.atoms[b]))
    xyz = np.vstack(
        [
            xyz,
            [
                [2, 2, 0],
                [2, 1.04, 0],
                [2.9, 2.2, 0],
                [20, 0, 0],
                [21, 0, 0],
                [30, 0, 0],
                [31, 1, 0],
            ],
        ]
    )
    mapping = RelayAcylAtomMap(acyl, 15, 16, 17, 18, 19, (20, 21))
    return topology, xyz, mapping


def build(relay: tuple, **kwargs: object) -> tuple:
    """Construct a synthetic bounded seed with explicit movement authorization."""
    topology, xyz, mapping = relay
    options = dict(
        active_atom_ids=tuple(range(1, len(xyz) + 1)),
        max_heavy_displacement_A=6.0,
        max_hydrogen_displacement_A=6.0,
        qm_charge=-1,
    )
    options.update(kwargs)
    return build_balanced_acyl_seed(topology, xyz, mapping, **options)


def test_balanced_local_product(relay: tuple) -> None:
    topology, source, mapping = relay
    original = source.copy()
    bonds = [(b.atom1.idx, b.atom2.idx) for b in topology.bonds]
    product, evidence = build(relay)
    np.testing.assert_array_equal(source, original)
    assert bonds == [(b.atom1.idx, b.atom2.idx) for b in topology.bonds]
    assert evidence["proton_inventory"] == dict(
        thr_n=2, thr_og=0, leaving_n=2, water_o=2, ash_oxygen=1
    )
    assert evidence["distances_A"] == pytest.approx(
        {"OG1_C": 1.4, "C_O": 1.23, "C_Nleaving": 3.0}
    )
    donor = evidence["donor_water_h_1based"] - 1
    assert donor == 15
    np.testing.assert_array_equal(product[4], source[donor])
    fixed = [0, 1, 2, 3, 8, 9, 14, 16, 17, 18, 19, 20]
    np.testing.assert_array_equal(product[fixed], source[fixed])
    assert np.linalg.norm(product[7] - product[11]) == pytest.approx(
        np.linalg.norm(source[7] - source[11])
    )
    assert np.linalg.norm(product[7] - product[donor]) == pytest.approx(1.03)
    assert evidence["qm_charge_before"] == evidence["qm_charge_after"] == -1
    assert len(evidence["source_coordinate_sha256"]) == 64
    second, report = build(relay)
    np.testing.assert_array_equal(second, product)
    assert report == evidence


@pytest.mark.parametrize(
    "change",
    [
        "shape",
        "nan",
        "cap",
        "cap_nan",
        "range",
        "charge",
        "outside",
        "water",
        "element",
        "asp_element",
        "ash_h",
        "asp_h",
        "thr_h",
        "duplicate",
        "index",
        "amide",
        "anchor",
        "proton_geometry",
    ],
)
def test_reject_invalid_chemistry_and_authorization(relay: tuple, change: str) -> None:
    topology, xyz, mapping = relay
    options = {}
    if change == "shape":
        xyz = xyz[:-1]
    elif change == "nan":
        xyz[0, 0] = np.nan
    elif change in ("cap", "cap_nan"):
        options["max_heavy_displacement_A"] = 0.01 if change == "cap" else np.nan
    elif change == "range":
        options["ester_distance"] = 4
    elif change == "charge":
        options["qm_charge"] = -1.5
    elif change == "outside":
        options["active_atom_ids"] = (1,)
    elif change == "water":
        topology.bonds.append(parmed.Bond(topology.atoms[14], topology.atoms[13]))
    elif change == "element":
        topology.atoms[14].atomic_number = 6
    elif change == "asp_element":
        topology.atoms[19].atomic_number = 6
    elif change == "ash_h":
        topology.bonds.append(parmed.Bond(topology.atoms[17], topology.atoms[12]))
    elif change == "asp_h":
        topology.bonds.append(parmed.Bond(topology.atoms[19], topology.atoms[12]))
    elif change == "thr_h":
        topology.bonds.append(parmed.Bond(topology.atoms[0], topology.atoms[12]))
    elif change == "duplicate":
        mapping = replace(mapping, water_h1=17)
    elif change == "index":
        mapping = replace(mapping, water_o=200)
    elif change == "amide":
        topology.atoms[11].atomic_number = 8
    elif change == "anchor":
        topology.atoms[10].atomic_number = 8
    elif change == "proton_geometry":
        xyz[4] = [9, 9, 9]
    with pytest.raises(ValueError):
        build((topology, xyz, mapping), **options)


def test_preserves_internal_components_with_distant_acyl_anchor(relay: tuple) -> None:
    topology, xyz, mapping = relay
    xyz[[5, 6, 10]] += [2, 0, 0]
    result, _ = build(relay)
    assert np.linalg.norm(result[5] - result[10]) == pytest.approx(
        np.linalg.norm(xyz[5] - xyz[10])
    )


def test_reports_close_contacts(relay: tuple) -> None:
    topology, xyz, mapping = relay
    product, _ = build(relay)
    xyz[13] = product[6] + [0.01, 0, 0]
    _, evidence = build(relay)
    assert evidence["close_contacts"]


def test_tip3p_hh_constraint_is_not_a_chemical_proton_parent(relay: tuple) -> None:
    topology, source, mapping = relay
    for index in (15, 16):
        topology.atoms[index].residue = topology.atoms[14].residue
    topology.bonds.append(parmed.Bond(topology.atoms[15], topology.atoms[16]))
    original = source.copy()
    bonds = [(b.atom1.idx, b.atom2.idx) for b in topology.bonds]
    result, evidence = build(relay)
    assert evidence["omitted_topology_constraints_1based"] == [[16, 17]]
    assert bonds == [(b.atom1.idx, b.atom2.idx) for b in topology.bonds]
    np.testing.assert_array_equal(source, original)
    np.testing.assert_array_equal(result[14], source[14])
    assert evidence["proton_inventory"]["water_o"] == 2


def test_reject_unrecognized_hydrogen_hydrogen_edge(relay: tuple) -> None:
    topology, _, _ = relay
    topology.bonds.append(parmed.Bond(topology.atoms[15], topology.atoms[16]))
    with pytest.raises(ValueError, match="Unrecognized water H-H"):
        build(relay)


def test_translate_only_and_bounded_search_preserve_fragments(relay: tuple) -> None:
    from nyl12_covalent.relay_seed_search import search_balanced_acyl_seed

    topology, source, mapping = relay
    original = source.copy()
    base, report = build(relay, idealize_ester_angle=False)
    np.testing.assert_allclose(base[10] - base[5], source[10] - source[5])
    assert report["idealize_ester_angle"] is False
    result, evidence = search_balanced_acyl_seed(
        topology,
        source,
        mapping,
        active_atom_ids=tuple(range(1, len(source) + 1)),
        qm_charge=-1,
    )
    assert evidence["candidate_count"] == 285
    assert evidence["water_oxygen_unchanged"]
    np.testing.assert_array_equal(source, original)
    np.testing.assert_array_equal(result[14], source[14])
    for a, b in ((5, 10), (7, 11), (7, 12)):
        assert np.linalg.norm(result[a] - result[b]) == pytest.approx(
            np.linalg.norm(source[a] - source[b])
        )
