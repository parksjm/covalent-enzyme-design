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
    stage2 = stage2_schedule(-1.7, 0.75, 1.3)
    assert len(stage2) == 21
    assert {row["attack_target_A"] for row in stage2} == {1.55}
    assert stage2[-1]["proton_difference_target_A"] == 0.8
    assert {row["relay_retention_target_A"] for row in stage2} == {0.75}
    assert {row["og1_hg1_minimum_A"] for row in stage2} == {1.3}


def test_relay_uses_od1_and_retains_the_first_transfer_coordinate() -> None:
    first = restraint_text(1, 1.55, 0.8, None)
    second = restraint_text(2, 1.55, 0.8, 0.75)
    assert first.count("&rst") == 2
    assert second.count("&rst") == 3
    assert "13916,13917,14559,13917" in first
    assert "14559,13917,13905,13917" in second
    assert "13916,13917,14559,13917" in second
    assert "r2=0.7500000000, r3=0.7500000000" in second
    assert "14560" not in first + second


@pytest.mark.parametrize(
    ("retention", "minimum", "message"),
    [
        (0.29, 1.3, "Relay-retention target"),
        (2.01, 1.3, "Relay-retention target"),
        (0.75, 1.19, "OG1-HG1 minimum"),
        (0.75, 2.01, "OG1-HG1 minimum"),
    ],
)
def test_stage2_rejects_invalid_proton_ownership_controls(
    retention: float, minimum: float, message: str
) -> None:
    """Reject invalid retained-coordinate or ownership thresholds."""
    with pytest.raises(ValueError, match=message):
        stage2_schedule(-1.7, retention, minimum)


@pytest.mark.parametrize(
    ("stage", "minimum", "message"),
    [
        (1, 0.75, "Stage I must not"),
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
