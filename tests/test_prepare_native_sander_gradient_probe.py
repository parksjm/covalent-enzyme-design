"""Tests for the one-cycle native-SANDER probe input preparation."""

import pytest

from scripts.prepare_native_sander_gradient_probe import prepare_input


def test_prepare_input_changes_only_temperature_mode_and_budget() -> None:
    source = """&qmmm
 qmmask='@1-3,5', qmcharge=-2, qm_theory='DFTB3', qmcut=10.0,
 scfconv=1.0d-10, dftb_telec=300, qmmm_int=1,
/
&dlfind
 dimer=.true., optalg='LBFGS', refcrd='/tmp/ref.ncrst',
 coords2(1)='/tmp/mode.ncrst', wtmask='@1-3', active=':1-3,5',
 tol=1.d-4, tole=1.d-5, maxcycle=200, maxstep=0.1,
/
"""

    prepared = prepare_input(source)

    assert "dftb_telec=0.0" in prepared
    assert "dimer=.false." in prepared
    assert "maxcycle=1" in prepared
    for unchanged in (
        "qmmask='@1-3,5'",
        "qmcharge=-2",
        "qm_theory='DFTB3'",
        "qmcut=10.0",
        "scfconv=1.0d-10",
        "active=':1-3,5'",
        "refcrd='/tmp/ref.ncrst'",
    ):
        assert unchanged in prepared


def test_prepare_input_rejects_unexpected_temperature() -> None:
    with pytest.raises(ValueError, match="dftb_telec"):
        prepare_input(
            "&qmmm qmmask='@1', qmcharge=0, qm_theory='DFTB3', qmcut=10, "
            "scfconv=1.d-10, dftb_telec=0 / &dlfind dimer=.true., "
            "optalg='LBFGS', active=':1', maxcycle=200 /"
        )
