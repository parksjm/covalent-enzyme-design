"""Geometric invariants and final movement gates of one acyl-arm torsion."""

from __future__ import annotations

import numpy as np
import pytest

from test_acylenzyme_seed import seed  # noqa: F401
from test_relay_acyl_seed import relay  # noqa: F401
from nyl12_covalent.relay_seed_refinement import refine_acyl_torsion_seed
from nyl12_covalent.relay_seed_search import search_balanced_acyl_seed


@pytest.fixture
def refinement(relay: tuple) -> tuple:  # noqa: F811
    topology, source, mapping = relay
    active = tuple(range(1, len(source) + 1))
    xyz, evidence = search_balanced_acyl_seed(
        topology, source, mapping, active_atom_ids=active, qm_charge=-1
    )
    return topology, source, mapping, xyz, evidence, active


def test_refinement_preserves_ester_leaving_water_and_internal_bonds(
    refinement: tuple,
) -> None:
    topology, source, mapping, base, evidence, active = refinement
    original, base_original = source.copy(), base.copy()
    result, report = refine_acyl_torsion_seed(
        topology,
        source,
        mapping,
        base,
        evidence,
        active_atom_ids=active,
        max_heavy_displacement_A=10,
        max_hydrogen_displacement_A=10,
    )
    summary = report["search"]["acyl_torsion"]
    assert summary["candidate_count"] == 24
    assert summary["torsion_axis_1based"] == [mapping.acyl.carbon, 11]
    assert summary["best_torsion_degrees"] in range(-180, 180, 15)
    np.testing.assert_array_equal(source, original)
    np.testing.assert_array_equal(base, base_original)
    np.testing.assert_array_equal(
        result[[3, 5, 6, 7, 10, 11, 12, 14, 16]], base[[3, 5, 6, 7, 10, 11, 12, 14, 16]]
    )
    assert report["proton_inventory"] == evidence["construction"]["proton_inventory"]
    assert (
        report["source_coordinate_sha256"]
        == evidence["construction"]["source_coordinate_sha256"]
    )
    assert np.linalg.norm(result[5] - result[10]) == pytest.approx(
        np.linalg.norm(base[5] - base[10])
    )
    assert report["distances_A"] == pytest.approx(
        evidence["construction"]["distances_A"]
        | {"C_Nleaving": evidence["best_choice"]["distance_A"]}
    )


@pytest.mark.parametrize(
    "bad", ["source_shape", "base_shape", "nan", "cap_zero", "cap_tight", "active"]
)
def test_refinement_rejects_invalid_geometry_or_caps(
    refinement: tuple, bad: str
) -> None:
    topology, source, mapping, base, evidence, active = refinement
    options = dict(
        active_atom_ids=active,
        max_heavy_displacement_A=10,
        max_hydrogen_displacement_A=10,
    )
    if bad == "source_shape":
        source = source[:-1]
    elif bad == "base_shape":
        base = base[:-1]
    elif bad == "nan":
        base[0, 0] = np.nan
    elif bad == "cap_zero":
        options["max_heavy_displacement_A"] = 0
    elif bad == "cap_tight":
        options["max_hydrogen_displacement_A"] = 0.001
    elif bad == "active":
        options["active_atom_ids"] = (1,)
    with pytest.raises(ValueError):
        refine_acyl_torsion_seed(topology, source, mapping, base, evidence, **options)
