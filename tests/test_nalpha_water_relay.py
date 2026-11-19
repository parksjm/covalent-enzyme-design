"""Check periodic relay geometry against explicit lattice reference searches."""

from itertools import product

import numpy as np
import pytest

from nyl12_covalent.nalpha_water_relay import (
    angle,
    candidate_rank,
    cell_matrix,
    close_waters,
    hbond,
    minimum_image,
    relay_class,
)


def test_minimum_image_skew_matches_exhaustive_lattice() -> None:
    """Fractional rounding alone must not classify octahedral-cell distances."""
    cell = cell_matrix(np.array([20.0, 20.0, 20.0]), np.array([109.471219] * 3))
    rng = np.random.default_rng(20261118)
    fractional = rng.uniform(-3, 3, (100, 3))
    vectors = fractional @ cell
    shifts = np.array(list(product(range(-5, 6), repeat=3))) @ cell
    expected = np.min(np.linalg.norm(vectors[:, None] - shifts, axis=2), axis=1)
    actual = np.linalg.norm(minimum_image(vectors, cell), axis=1)
    np.testing.assert_allclose(actual, expected, atol=1e-10)
    rounded = (fractional - np.rint(fractional)) @ cell
    assert np.any(np.linalg.norm(rounded, axis=1) > expected + 1e-5)
    assert minimum_image(np.empty((0, 3)), cell).shape == (0, 3)
    with pytest.raises(ValueError, match="finite"):
        minimum_image(np.array([[np.nan, 0, 0]]), cell)


def test_near_water_filter_matches_exhaustive_distances() -> None:
    """Conservative filtering retains every actual periodic water contact."""
    cell = cell_matrix(np.array([20.0, 20.0, 20.0]), np.array([109.471219] * 3))
    positions = np.random.default_rng(10).uniform(-2, 2, (1000, 3)) @ cell
    anchor = np.array([3.0, 4.0, 5.0])
    expected = np.flatnonzero(
        np.linalg.norm(minimum_image(positions - anchor, cell), axis=1) <= 3.5
    )
    np.testing.assert_array_equal(close_waters(positions, anchor, cell, 3.5), expected)
    with pytest.raises(ValueError, match="Cell too small"):
        close_waters(positions, anchor, np.eye(3), 3.5)


def test_hbond_periodic_direction_and_thresholds() -> None:
    """A close water with a hydrogen pointing away cannot pass the relay leg."""
    cell = np.eye(3) * 20
    xyz = np.array(
        [[19.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.8, 0.0, 0.0], [18.0, 0.0, 0.0]]
    )
    passing = hbond(xyz, 0, 1, 2, cell)
    assert passing["pass"]
    assert passing["ha_A"] == pytest.approx(1.8)
    assert passing["angle_deg"] == pytest.approx(180)
    assert not hbond(xyz, 0, 3, 2, cell)["pass"]
    with pytest.raises(ValueError, match="Degenerate"):
        angle(np.zeros(3), np.zeros(3), np.ones(3))


def test_classification_does_not_reject_missing_future_proton() -> None:
    """Keep leaving-leg-ready neutral NH2 frames distinct from strict two-leg."""
    assert relay_class(True, True) == "strict_both_leg_hbond"
    assert relay_class(False, True) == "leaving_leg_hbond"
    assert relay_class(False, False) == "heavy_bridge_only"
    assert relay_class(True, False) == "heavy_bridge_only"
    common = {
        "tier": "A",
        "water_h_to_nleave_A": 2.0,
        "water_h_to_nleave_angle_deg": 160.0,
        "native_thr_rmsd_A": 0.5,
        "case": "r1",
        "frame": 1,
        "water_resid": 1500,
    }
    strict = {**common, "relay_class": "strict_both_leg_hbond"}
    leaving = {**common, "relay_class": "leaving_leg_hbond"}
    heavy = {**common, "relay_class": "heavy_bridge_only"}
    assert candidate_rank(strict) < candidate_rank(leaving) < candidate_rank(heavy)


def test_native_thr_gate_preserves_rigid_alignment_and_rejects_og1_shift() -> None:
    """Native gate must ignore global pose and detect held-out OG1 displacement."""
    from pathlib import Path
    from types import SimpleNamespace

    from nyl12_covalent.catalytic_alignment import SITE_RESIDUES, _cif_atoms
    from nyl12_covalent.nalpha_water_relay import native_geometry

    reference = _cif_atoms(
        Path(__file__).resolve().parents[1] / "cif/n12c.cif",
        "A",
        {auth for _, auth, _ in SITE_RESIDUES},
    )
    residues = [None] * 970
    points = []
    mapped = {}
    for amber, auth, name in SITE_RESIDUES:
        atoms = []
        for (number, atom_name), (_, point) in reference.items():
            if number != auth:
                continue
            index = len(points)
            atoms.append(SimpleNamespace(idx=index, name=atom_name, atomic_number=6))
            points.append(point)
            mapped[(auth, atom_name)] = index
        residues[amber - 1] = SimpleNamespace(
            name="ASH" if amber == 970 else name, atoms=atoms
        )
    topology = SimpleNamespace(residues=residues)
    xyz = (
        np.array(points)
        @ np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        + 15
    )
    cell = np.eye(3) * 200
    native = native_geometry(xyz, topology, reference, cell, mapped[(233, "N")])
    assert native["native_pass"]
    assert native["native_thr_rmsd_A"] == pytest.approx(0, abs=1e-10)
    xyz[mapped[(233, "OG1")]] += np.array([2.0, 0.0, 0.0])
    assert not native_geometry(xyz, topology, reference, cell, mapped[(233, "N")])[
        "native_pass"
    ]
    residues[969].name = "LYS"
    with pytest.raises(ValueError, match="mapping mismatch"):
        native_geometry(xyz, topology, reference, cell, mapped[(233, "N")])
