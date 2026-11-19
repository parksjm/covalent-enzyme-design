"""Distinguish a TS-like shared transfer H from spectator/acid/overlap defects."""

from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.neb_spacing_audit import spacing_image_audit


def model_and_source() -> tuple[SimpleNamespace, np.ndarray]:
    """Return a small explicitly labeled donor/acceptor and spectator inventory."""
    xyz = np.array(
        [
            [0.0, 0, 0],
            [3, 0, 0],
            [10, 0, 0],
            [20, 0, 0],
            [1, 0, 0],
            [11, 0, 0],
            [21, 0, 0],
            [31, 0, 0],
            [30, 0, 0],
        ]
    )
    atoms = [SimpleNamespace(atomic_number=z) for z in (8, 7, 8, 8, 1, 1, 1, 1, 7)]
    model = SimpleNamespace(
        qm=np.arange(9),
        hydrogens=np.array([4, 5, 6, 7]),
        topology=SimpleNamespace(atoms=atoms),
        atom_map=dict(water_h1=5, thr_og=1, leaving_n=2, ash_proton=7, ash_oxygen=4),
    )
    return model, xyz


def test_normal_and_transferred_H_preserve_spectators() -> None:
    model, source = model_and_source()
    assert spacing_image_audit(model, source, source)["strict_H_gate_passed"]
    product = source.copy()
    product[4] = [2, 0, 0]
    report = spacing_image_audit(model, product, source)
    assert report["strict_H_gate_passed"]
    assert report["spectator_parent_changes"] == []


def test_shared_transfer_is_explicitly_flagged_and_not_automatically_accepted() -> None:
    model, source = model_and_source()
    candidate = source.copy()
    candidate[4] = [1.5, 0, 0]
    report = spacing_image_audit(model, candidate, source)
    assert report["strict_H_gate_passed"] is False
    assert report["shared_transfer20318_only_requires_root_decision"] is True
    assert report["no_compressed_H_contacts"] is True
    offender = report["offending_hydrogens"][0]
    assert offender["hydrogen_1based"] == 5
    assert {offender["nearest_heavy_1based"], offender["second_heavy_1based"]} == {1, 2}
    assert offender["nearest_distance_A"] == pytest.approx(1.5)


@pytest.mark.parametrize("defect", ["spectator", "acid", "overlap", "unassigned"])
def test_nonreactive_defects_are_never_shared_transfer_exception(defect: str) -> None:
    model, source = model_and_source()
    candidate = source.copy()
    if defect == "spectator":
        candidate[7] = [19, 0, 0]
    elif defect == "acid":
        model.atom_map["ash_oxygen"] = 2
    elif defect == "overlap":
        candidate[4] = [0.1, 0, 0]
    elif defect == "unassigned":
        candidate[7] = [40, 0, 0]
    report = spacing_image_audit(model, candidate, source)
    assert report["strict_H_gate_passed"] is False
    assert report["shared_transfer20318_only_requires_root_decision"] is False
