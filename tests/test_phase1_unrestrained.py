"""Check that releasing phase-I biases does not replace the QM Hamiltonian."""

import pytest

from provenance.phase1_unrestrained_20261025.run_optimization import (
    unrestrained_input,
)


SOURCE = """phase-I
&cntrl
 nmropt=1, ntr=1, restraint_wt=10.0, restraintmask='!(@1-98)',
 ntc=1, ntf=1, ntb=1, ntwx=5,
/
&qmmm
 qm_theory='DFTB3', qmcharge=-2, qmmask='@1-98',
/
&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint_values.dat
&dlfind optalg='LBFGS', active='@*', tol=1.d-4, maxcycle=1000, /
"""


def test_release_preserves_qm_and_has_no_frozen_atoms() -> None:
    """Only artificial biases and the declared cycle budget may change."""
    result = unrestrained_input(SOURCE)
    assert (
        result.split("&qmmm")[1].split("/")[0]
        == (SOURCE.split("&qmmm")[1].split("/")[0])
    )
    assert "ntr=0" in result and "nmropt=0" in result
    assert "ntc=1, ntf=1, ntb=1, ntwx=5" in result
    assert "restraint_wt" not in result and "restraintmask" not in result
    assert "DISANG" not in result and "DUMPAVE" not in result
    assert "active='@*'" in result and "tol=1.d-4" in result
    assert "maxcycle=2000" in result


@pytest.mark.parametrize(
    "original,replacement",
    [
        ("nmropt=1", "nmropt=0"),
        ("ntr=1", "ntr=0"),
        ("restraint_wt=10.0", "restraint_wt=5.0"),
        ("active='@*'", "active='@1-98'"),
        ("maxcycle=1000", "maxcycle=10"),
        ("&wt type='END' /", ""),
    ],
)
def test_unexpected_input_fails(original: str, replacement: str) -> None:
    """Reject controls that were not scientifically reviewed."""
    with pytest.raises(ValueError):
        unrestrained_input(SOURCE.replace(original, replacement))
