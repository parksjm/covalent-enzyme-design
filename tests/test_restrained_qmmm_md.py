"""Tests for the restrained DFTB3/MM MD pilot inputs."""

from nyl12_covalent.restrained_qmmm_md import (
    heating_input,
    minimization_input,
    proton_restraint,
)

AUDITED = """accepted
&qmmm
 qmmask='@1-3', qmcharge=-2, qm_theory='DFTB3',
/
"""


def test_proton_restraint_is_flat_bottom_difference_coordinate() -> None:
    text = proton_restraint()
    assert "iat=13916,13917,13905,13917" in text
    assert "rstwt=1,-1,0,0" in text
    assert "r2=0.600000" in text and "r3=0.850000" in text
    assert "rk2=50.0, rk3=50.0" in text


def test_minimization_has_no_positional_restraints() -> None:
    text = minimization_input(AUDITED)
    assert "imin=1" in text and "ntr=0" in text and "nmropt=1" in text
    assert "restraint_wt" not in text and "restraintmask" not in text
    assert "qm_theory='DFTB3'" in text


def test_heating_is_one_ps_nvt_with_flexible_hydrogens() -> None:
    text = heating_input(AUDITED)
    assert "nstlim=2000, dt=0.0005" in text
    assert "ntb=1, ntp=0" in text
    assert "ntc=1, ntf=1" in text
    assert "value1=25.0, value2=300.0" in text
