import pytest

from nyl12_covalent.unrestrained_relay import (
    driven_input,
    es_input,
    restraint_text,
    schedule,
    stage2_schedule,
)


BASE = """Title
&cntrl
 nmropt=0, imin=1, ntr=0, ntwx=5,
/
&qmmm
 qm_theory='DFTB3',
/
&dlfind optalg='LBFGS', maxcycle=2000, verbosity=2, /
"""


def test_inputs_have_no_positional_restraints() -> None:
    relaxed = es_input(BASE)
    driven = driven_input(BASE)
    for text in (relaxed, driven):
        assert "ntr=0" in text
        assert "restraint_wt" not in text
        assert "restraintmask" not in text
        assert "verbosity=1" in text
    assert "nmropt=1" in driven
    assert "DISANG=restraints.RST" in driven
    assert driven.index("&wt type='END'") < driven.index("&dlfind")


def test_two_stage_schedule_does_not_extrapolate_bonds() -> None:
    initial = {
        "OG1_C": 3.2,
        "OG1_HG1": 1.0,
        "Asp272_OD1_HG1": 1.8,
    }
    stage1 = schedule(initial)
    assert len(stage1) == 21
    assert stage1[0]["attack_target_A"] == 3.2
    assert stage1[-1]["attack_target_A"] == 1.55
    assert stage1[-1]["proton_difference_target_A"] == 0.8
    stage2 = stage2_schedule(-1.7, 1.35)
    assert len(stage2) == 21
    assert {row["attack_target_A"] for row in stage2} == {1.55}
    assert stage2[-1]["proton_difference_target_A"] == 0.8
    assert {row["donor_separation_minimum_A"] for row in stage2} == {1.35}


def test_relay_uses_od1_and_explicit_donor_separation_wall() -> None:
    first = restraint_text(1, 1.55, 0.8, None)
    second = restraint_text(2, 1.55, 0.8, 1.35)
    assert first.count("&rst") == 2
    assert second.count("&rst") == 3
    assert "13916,13917,14559,13917" in first
    assert "14559,13917,13905,13917" in second
    assert "13916,13917,0,0,0,0,0,0" in second
    assert "r2=1.3500000000, r3=10.0000000000" in second
    assert "14560" not in first + second


@pytest.mark.parametrize("minimum", [1.19, 2.01])
def test_stage2_rejects_invalid_donor_separation_wall(minimum: float) -> None:
    """Reject a wall outside the declared chemically useful interval."""
    with pytest.raises(ValueError, match="between 1.2 and 2.0"):
        stage2_schedule(-1.7, minimum)


@pytest.mark.parametrize(
    ("stage", "minimum", "message"),
    [
        (1, 1.35, "Stage I must not"),
        (2, None, "Stage II requires"),
        (3, None, "stage must be 1 or 2"),
    ],
)
def test_restraints_require_stage_specific_donor_wall(
    stage: int, minimum: float | None, message: str
) -> None:
    """Fail rather than silently omit or misapply the donor wall."""
    with pytest.raises(ValueError, match=message):
        restraint_text(stage, 1.55, 0.8, minimum)
