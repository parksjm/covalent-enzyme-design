"""Chemical invariants of atom-preserving acylenzyme coordinate hypotheses."""

from dataclasses import replace

import numpy as np
import parmed
import pytest

from nyl12_covalent.acylenzyme_seed import AcylAtomMap, build_acylenzyme_seed


@pytest.fixture
def seed() -> tuple[parmed.Structure, np.ndarray, AcylAtomMap]:
    structure = parmed.Structure()
    names = [
        "N",
        "H1",
        "H2",
        "OG1",
        "HG1",
        "C",
        "O",
        "N",
        "CB",
        "CA",
        "C5",
        "C1",
        "HN",
        "O",
    ]
    elements = [7, 1, 1, 8, 1, 6, 8, 7, 6, 6, 6, 6, 1, 8]
    for i, (name, element) in enumerate(zip(names, elements, strict=True)):
        residue = "THR" if i in (0, 1, 2, 3, 4, 8, 9) else "P6I"
        # ParmEd residues must be contiguous, so explicitly reuse the Thr object.
        structure.add_atom(
            parmed.Atom(name=name, atomic_number=element), residue, i + 1
        )
    thr = structure.atoms[0].residue
    for index in (1, 2, 3, 4, 8, 9):
        structure.atoms[index].residue = thr
    bonds = [
        (0, 1),
        (0, 2),
        (0, 9),
        (9, 8),
        (8, 3),
        (3, 4),
        (5, 6),
        (5, 7),
        (5, 10),
        (7, 11),
        (7, 12),
    ]
    for a, b in bonds:
        structure.bonds.append(parmed.Bond(structure.atoms[a], structure.atoms[b]))
    coordinates = np.array(
        [
            [-1, 0, 0],
            [-2, 0, 0],
            [-1, 1, 0],
            [0, 0, 0],
            [-1, 0, 1],
            [1.6, 0, 0],
            [1.8, 1.3, 0],
            [2.2, -1.2, 0],
            [0, 0, -1.4],
            [-1, 0, -1.4],
            [2, 0.2, 1.4],
            [3, -1.8, 0.5],
            [1.7, -1.8, -0.7],
            [30, 30, 30],
        ],
        dtype=float,
    )
    return structure, coordinates, AcylAtomMap(1, 2, 3, 4, 5, 6, 7, 8)


def test_planar_product_preserves_identity_and_rigid_fragments(seed: tuple) -> None:
    topology, source, mapping = seed
    original = source.copy()
    bonds = [(b.atom1.idx, b.atom2.idx) for b in topology.bonds]
    names = [a.name for a in topology.atoms]
    product, report = build_acylenzyme_seed(topology, source, mapping)
    np.testing.assert_array_equal(source, original)
    assert names == [a.name for a in topology.atoms]
    assert bonds == [(b.atom1.idx, b.atom2.idx) for b in topology.bonds]
    assert report["distances_A"]["OG1_C"] == pytest.approx(1.4)
    assert report["distances_A"]["C_O"] == pytest.approx(1.23)
    assert report["distances_A"]["C_Nleaving"] == pytest.approx(2.85)
    assert report["distances_A"]["Nleaving_HG1"] == pytest.approx(1.03)
    np.testing.assert_array_equal(
        product[[0, 1, 2, 3, 8, 9, 13]], source[[0, 1, 2, 3, 8, 9, 13]]
    )
    for a, b in ((5, 10), (7, 11), (7, 12), (11, 12)):
        assert np.linalg.norm(product[a] - product[b]) == pytest.approx(
            np.linalg.norm(source[a] - source[b])
        )
    vectors = product[[3, 6, 10]] - product[5]
    assert abs(np.linalg.det(vectors)) < 1e-12
    directions = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    assert np.dot(directions[0], directions[1]) == pytest.approx(
        np.dot(directions[1], directions[2])
    )
    assert 95 <= report["OG1_C_acyl_angle_degrees"] <= 140
    np.testing.assert_array_equal(product[10], source[10])
    assert set(report["moved_indices_1based"]) == {5, 6, 7, 8, 12, 13}


def test_reports_nonbonded_clash(seed: tuple) -> None:
    topology, source, mapping = seed
    product, _ = build_acylenzyme_seed(topology, source, mapping)
    source[13] = product[6] + [0.1, 0, 0]
    _, report = build_acylenzyme_seed(topology, source, mapping)
    assert any(
        contact["atoms_1based"] == [7, 14] for contact in report["close_contacts"]
    )


def test_explicit_preparation_distances(seed: tuple) -> None:
    topology, source, mapping = seed
    _, report = build_acylenzyme_seed(
        topology, source, mapping, ester_distance=1.65, leaving_distance=2.3
    )
    assert report["distances_A"]["OG1_C"] == pytest.approx(1.65)
    assert report["distances_A"]["C_Nleaving"] == pytest.approx(2.3)
    with pytest.raises(ValueError, match="distances"):
        build_acylenzyme_seed(topology, source, mapping, leaving_distance=1.5)


def test_rigid_orientation_reduces_overlap_reproducibly(seed: tuple) -> None:
    topology, source, mapping = seed
    plain, _ = build_acylenzyme_seed(topology, source, mapping)
    source[13] = plain[11] + [0.05, 0, 0]
    oriented, report = build_acylenzyme_seed(
        topology, source, mapping, orient_leaving=True
    )
    repeated, _ = build_acylenzyme_seed(topology, source, mapping, orient_leaving=True)
    np.testing.assert_array_equal(oriented, repeated)
    assert np.linalg.norm(oriented[11] - source[13]) > 0.05
    assert abs(report["leaving_rigid_rotation"]["degrees"]) <= 15
    assert np.linalg.norm(oriented[7] - oriented[4]) == pytest.approx(1.03)
    for a, b in ((7, 11), (7, 12), (11, 12)):
        assert np.linalg.norm(oriented[a] - oriented[b]) == pytest.approx(
            np.linalg.norm(source[a] - source[b])
        )


@pytest.mark.parametrize(
    "change",
    [
        "names",
        "elements",
        "index",
        "duplicate",
        "bond",
        "cycle",
        "proton",
        "nan",
        "shape",
        "substituent",
        "degenerate",
        "residue",
    ],
)
def test_rejects_invalid_source(seed: tuple, change: str) -> None:
    topology, source, mapping = seed
    if change == "names":
        topology.atoms[3].name = "OX"
    elif change == "elements":
        topology.atoms[3].atomic_number = 6
    elif change == "index":
        mapping = replace(mapping, carbon=100)
    elif change == "duplicate":
        mapping = replace(mapping, carbon=mapping.oxygen)
    elif change == "bond":
        topology.bonds[-1].delete()
        topology.bonds.append(parmed.Bond(topology.atoms[4], topology.atoms[7]))
    elif change == "cycle":
        topology.bonds.append(parmed.Bond(topology.atoms[10], topology.atoms[11]))
    elif change == "proton":
        source[4] = [10, 10, 10]
    elif change == "nan":
        source[0, 0] = np.nan
    elif change == "shape":
        source = source[:-1]
    elif change == "substituent":
        topology.atoms[10].atomic_number = 8
    elif change == "degenerate":
        source[10] = source[5]
    elif change == "residue":
        topology.atoms[0].residue.name = "SER"
    with pytest.raises(ValueError):
        build_acylenzyme_seed(topology, source, mapping)


def test_explicit_leaving_n_source_state(seed: tuple) -> None:
    topology, source, mapping = seed
    source[4] = source[7] + [1.03, 0, 0]
    original = source.copy()
    _, report = build_acylenzyme_seed(
        topology, source, mapping, source_proton_state="leaving_n"
    )
    assert report["source_proton_state"] == "leaving_n"
    assert report["distances_A"]["Nleaving_HG1"] == pytest.approx(1.03)
    np.testing.assert_array_equal(source, original)
    with pytest.raises(ValueError, match="proton state"):
        build_acylenzyme_seed(topology, source, mapping)


@pytest.mark.parametrize("state", ["leaving_n", "unknown"])
def test_rejects_mismatched_explicit_state(seed: tuple, state: str) -> None:
    topology, source, mapping = seed
    with pytest.raises(ValueError, match="proton state"):
        build_acylenzyme_seed(topology, source, mapping, source_proton_state=state)


@pytest.mark.parametrize("bend", [-30, 20])
def test_explicit_leaving_bend_preserves_product_distances(
    seed: tuple, bend: float
) -> None:
    topology, source, mapping = seed
    product, report = build_acylenzyme_seed(
        topology, source, mapping, leaving_bend_degrees=bend
    )
    assert report["leaving_bend_degrees"] == bend
    assert report["distances_A"]["C_Nleaving"] == pytest.approx(2.85)
    assert np.linalg.norm(product[7] - product[11]) == pytest.approx(
        np.linalg.norm(source[7] - source[11])
    )


@pytest.mark.parametrize("bend", [31, float("nan")])
def test_invalid_leaving_bend_rejected(seed: tuple, bend: float) -> None:
    topology, source, mapping = seed
    with pytest.raises(ValueError, match="Leaving bend"):
        build_acylenzyme_seed(topology, source, mapping, leaving_bend_degrees=bend)
