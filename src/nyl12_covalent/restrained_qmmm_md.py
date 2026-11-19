"""Inputs for restrained DFTB3/MM environmental relaxation."""

from __future__ import annotations

import re

N_ATOM = 13905
OG1_ATOM = 13916
HG1_ATOM = 13917


def qmmm_block(audited_input: str) -> str:
    """Extract the reviewed QM/MM Hamiltonian block from an accepted input."""
    match = re.search(r"&qmmm\b.*?\n/", audited_input, flags=re.DOTALL | re.IGNORECASE)
    if match is None:
        raise ValueError("Accepted input lacks one complete qmmm block")
    block = match.group(0)
    if "qm_theory='DFTB3'" not in block or "qmcharge=-2" not in block:
        raise ValueError("Expected the reviewed charge -2 DFTB3 Hamiltonian")
    return block


def proton_restraint(lower: float = 0.60, upper: float = 0.85) -> str:
    """Return a flat-bottom transferred-proton difference-coordinate restraint."""
    if not -1.0 < lower < upper < 2.0:
        raise ValueError("Invalid proton-coordinate flat-bottom bounds")
    return (
        f"&rst iat={OG1_ATOM},{HG1_ATOM},{N_ATOM},{HG1_ATOM},0,0,0,0, "
        f"rstwt=1,-1,0,0, r1={lower - 1.0:.6f}, r2={lower:.6f}, "
        f"r3={upper:.6f}, r4={upper + 1.0:.6f}, "
        "rk2=50.0, rk3=50.0, /\n"
    )


def minimization_input(audited_input: str, cycles: int = 500) -> str:
    """Build a short all-atom restrained QM/MM minimization input."""
    if cycles < 1:
        raise ValueError("cycles must be positive")
    return f"""Transferred-proton restrained all-atom DFTB3/MM pre-MD minimization
&cntrl
 imin=1, maxcyc={cycles}, ncyc=100, ntmin=1, dx0=0.01, drms=0.05,
 ntb=1, igb=0, ntc=1, ntf=1, cut=10.0, ifqnt=1, jfastw=4,
 ntr=0, nmropt=1, ntpr=25, ntwx=25, ntwr=100, ioutfm=1, ntxo=2,
/
{qmmm_block(audited_input)}
&wt type='END' /
DISANG=proton_coordinate.RST
DUMPAVE=proton_coordinate.dat
"""


def heating_input(audited_input: str, steps: int = 2_000) -> str:
    """Build a 25-to-300 K NVT restrained QM/MM heating pilot."""
    if steps < 1:
        raise ValueError("steps must be positive")
    return f"""Transferred-proton restrained DFTB3/MM 25-to-300 K NVT pilot
&cntrl
 imin=0, irest=0, ntx=1, nstlim={steps}, dt=0.0005,
 ntb=1, ntp=0, igb=0, cut=10.0, iwrap=1,
 ntt=3, gamma_ln=2.0, tempi=25.0, temp0=300.0, ig=-1,
 ntc=1, ntf=1, ifqnt=1, jfastw=4, ntr=0, nmropt=1,
 ntpr=100, ntwx=100, ntwr=500, ioutfm=1, ntxo=2,
/
{qmmm_block(audited_input)}
&wt type='TEMP0', istep1=0, istep2={steps}, value1=25.0, value2=300.0 /
&wt type='END' /
DISANG=proton_coordinate.RST
DUMPAVE=proton_coordinate.dat
"""
