"""Tests for the direct Thr233 proton-transfer construction."""

import pytest

from nyl12_covalent.direct_proton_transfer import (
    driven_input,
    release_input,
    restraint_text,
    schedule,
)

BASE = """Accepted 12 A input
&cntrl
 nmropt=0, ntr=0, maxcyc=10000,
/
&qmmm
 qm_theory='DFTB3', qmmask='@1-3',
/
&dlfind active='(:1)<:12.0', maxcycle=10000, /
"""


def test_schedule_is_bounded_and_uses_small_increments() -> None:
    rows = schedule(-1.09, final_difference=0.75, maximum_increment=0.10)
    targets = [row["target_A"] for row in rows]
    assert targets[-1] == pytest.approx(0.75)
    assert targets[0] > -1.09
    assert max(b - a for a, b in zip([-1.09, *targets[:-1]], targets)) <= 0.10


def test_driven_input_has_only_reaction_coordinate_restraint() -> None:
    text = driven_input(BASE)
    assert "nmropt=1" in text
    assert "ntr=0" in text
    assert "DISANG=restraints.RST" in text
    assert "active='(:1)<:12.0'" in text
    assert "restraint_wt" not in text
    assert "restraintmask" not in text
    assert "maxcycle=2000" in text


def test_release_input_removes_all_reaction_bias() -> None:
    text = release_input(BASE)
    assert "nmropt=0" in text
    assert "DISANG" not in text
    assert "ntr=0" in text
    assert "active='(:1)<:12.0'" in text


def test_restraint_is_direct_og1_to_n_difference_coordinate() -> None:
    text = restraint_text(0.75)
    assert "iat=13916,13917,13905,13917" in text
    assert "rstwt=1,-1,0,0" in text
    assert "r2=0.7500000000" in text
    assert text.count("&rst") == 1
