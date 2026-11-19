"""DFTB3/MM proton-transfer umbrella-sampling input builders."""

from __future__ import annotations

from nyl12_covalent.restrained_qmmm_md import HG1_ATOM, N_ATOM, OG1_ATOM, qmmm_block

CARBONYL_C_ATOM = 19916


def umbrella_restraint(target: float, force_constant: float = 50.0) -> str:
    """Return a harmonic difference-coordinate umbrella restraint."""
    if not -2.0 < target < 2.0 or force_constant <= 0.0:
        raise ValueError("Invalid umbrella target or force constant")
    return (
        f"&rst iat={OG1_ATOM},{HG1_ATOM},{N_ATOM},{HG1_ATOM},0,0,0,0, "
        f"rstwt=1,-1,0,0, r1={target - 1.0:.8f}, r2={target:.8f}, "
        f"r3={target:.8f}, r4={target + 1.0:.8f}, "
        f"rk2={force_constant:.1f}, rk3={force_constant:.1f}, /\n"
    )


def asynchronous_umbrella_restraint(target: float, force_constant: float = 50.0) -> str:
    """Return an umbrella on proton transfer minus the OG1 attack distance."""
    if not -6.0 < target < 2.0 or force_constant <= 0.0:
        raise ValueError("Invalid asynchronous umbrella target or force constant")
    return (
        f"&rst iat={OG1_ATOM},{HG1_ATOM},{N_ATOM},{HG1_ATOM},"
        f"{OG1_ATOM},{CARBONYL_C_ATOM},0,0, rstwt=1,-1,-1,0, "
        f"r1={target - 1.0:.8f}, r2={target:.8f}, "
        f"r3={target:.8f}, r4={target + 1.0:.8f}, "
        f"rk2={force_constant:.1f}, rk3={force_constant:.1f}, /\n"
    )


def md_input(
    audited_input: str,
    stage: str,
    steps: int,
    time_step_ps: float = 0.0005,
) -> str:
    """Build one NVT umbrella heating, equilibration, or production input."""
    if stage not in {"heating", "equilibration", "production"}:
        raise ValueError("Unsupported umbrella stage")
    if steps < 1 or time_step_ps <= 0.0:
        raise ValueError("Positive steps and time step are required")
    starting = (
        "irest=0, ntx=1, tempi=25.0"
        if stage == "heating"
        else "irest=1, ntx=5, tempi=300.0"
    )
    trajectory_stride = 500 if stage != "production" else 100
    weights = (
        f"&wt type='TEMP0', istep1=0, istep2={steps}, value1=25.0, value2=300.0 /\n"
        if stage == "heating"
        else ""
    )
    return f"""DFTB3/MM umbrella {stage}
&cntrl
 imin=0, {starting}, nstlim={steps}, dt={time_step_ps:.7f},
 ntb=1, ntp=0, igb=0, cut=10.0, iwrap=1,
 ntt=3, gamma_ln=2.0, temp0=300.0, ig=-1,
 ntc=1, ntf=1, ifqnt=1, jfastw=4, ntr=0, nmropt=1,
 ntpr=100, ntwx={trajectory_stride}, ntwr=500, ioutfm=1, ntxo=2,
/
{qmmm_block(audited_input)}
{weights}&wt type='DUMPFREQ', istep1=10 /
&wt type='END' /
DISANG=umbrella.RST
DUMPAVE=umbrella_coordinate.dat
"""
