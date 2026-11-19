import pytest

from nyl12_covalent.es_classical_relaxation import (
    minimization_input,
    solvent_relaxation_input,
)


def test_minimization_is_fully_unrestrained_and_saves_trajectory() -> None:
    text = minimization_input()
    assert "ntr=0" in text
    assert "nmropt=0" in text
    assert "ntwx=25" in text
    assert "restraint_wt" not in text
    assert "restraintmask" not in text


def test_minimization_rejects_nonpositive_cycles() -> None:
    with pytest.raises(ValueError, match="positive"):
        minimization_input(0)


def test_minimization_uses_explicit_case_label() -> None:
    assert minimization_input(case_label="WT frame 1553").startswith("WT frame 1553")
    with pytest.raises(ValueError, match="case_label"):
        minimization_input(case_label=" ")


def test_minimization_uses_requested_cycle_budget() -> None:
    text = minimization_input(max_cycles=20_000)
    assert "maxcyc=20000" in text
    assert "ncyc=10000" in text


def test_solvent_relaxation_restraints_are_explicit_and_scoped() -> None:
    text = solvent_relaxation_input(case_label="Template Model 7")
    assert text.startswith("Template Model 7")
    assert "ntr=1" in text
    assert "restraintmask=':1-1313 & !@H='" in text
    assert "restraint_wt=10.0" in text
    assert "nmropt=0" in text


def test_solvent_relaxation_rejects_empty_label() -> None:
    with pytest.raises(ValueError, match="case_label"):
        solvent_relaxation_input(case_label=" ")
