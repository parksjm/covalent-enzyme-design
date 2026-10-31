"""Amber MD protocol text for modular PA6/Nyl12 simulations."""

from __future__ import annotations

from pathlib import Path


STAGES = (
    "01_min_restrained",
    "02_min_unrestrained",
    "03_heat_nvt",
    "04_density_npt",
    "05_equil_weak_npt",
    "06_equil_npt",
    "07_prod_100ns",
)


def mdin_files(seed: int) -> dict[str, str]:
    """Return the seven Amber input files for one 100 ns replica.

    Args:
        seed: Distinct positive Amber random seed for heating.

    Returns:
        Mapping from input filename to complete mdin text.
    """
    if seed <= 0:
        raise ValueError("Amber seed must be positive")
    restraint = "!:WAT,Na+,Cl- & !@H="
    return {
        "01_min_restrained.in": _minimization(restraint),
        "02_min_unrestrained.in": _minimization(None),
        "03_heat_nvt.in": _heating(seed, restraint),
        "04_density_npt.in": _equilibration("density", 5.0, restraint),
        "05_equil_weak_npt.in": _equilibration("weak", 1.0, restraint),
        "06_equil_npt.in": _equilibration("unrestrained", None, None),
        "07_prod_100ns.in": _production(seed + 10_000),
    }


def replica_runner(pmemd_cuda: Path) -> str:
    """Return a restart-aware shell runner for one replica.

    Args:
        pmemd_cuda: Absolute path to the pinned PMEMD CUDA executable.

    Returns:
        Bash script that executes all unfinished stages in order.
    """
    if not pmemd_cuda.is_absolute():
        raise ValueError("pmemd.cuda path must be absolute")
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f'PMEMD_CUDA="{pmemd_cuda}"',
        'SYSTEM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"',
        'PRMTOP="${SYSTEM_DIR}/md_solvated.prmtop"',
        'INPCRD="${SYSTEM_DIR}/md_solvated.inpcrd"',
        "",
    ]
    stage_inputs = (
        ("01_min_restrained", "${INPCRD}", "${INPCRD}"),
        ("02_min_unrestrained", "01_min_restrained.rst7", None),
        ("03_heat_nvt", "02_min_unrestrained.rst7", "02_min_unrestrained.rst7"),
        ("04_density_npt", "03_heat_nvt.rst7", "03_heat_nvt.rst7"),
        ("05_equil_weak_npt", "04_density_npt.rst7", "04_density_npt.rst7"),
        ("06_equil_npt", "05_equil_weak_npt.rst7", None),
        ("07_prod_100ns", "06_equil_npt.rst7", None),
    )
    for stage, coordinates, reference in stage_inputs:
        lines.extend(
            (
                f'if [[ ! -s "{stage}.rst7" ]]; then',
                f'  "${{PMEMD_CUDA}}" -O -i {stage}.in -o {stage}.out \\',
                f'    -p "${{PRMTOP}}" -c "{coordinates}" -r {stage}.rst7 \\',
            )
        )
        final = f"    -inf {stage}.mdinfo"
        if stage.startswith(("03", "04", "05", "06", "07")):
            final += f" -x {stage}.nc"
        if reference is not None:
            final += f' -ref "{reference}"'
        lines.extend((final, "fi", ""))
    return "\n".join(lines)


def cpptraj_analysis() -> str:
    """Return modular-residue trajectory analysis commands."""
    return """parm ../md_solvated.prmtop
trajin 07_prod_100ns.nc
autoimage
rms first :1-1312@CA,C,N,O
atomicfluct out substrate_rmsf_by_atom.dat :1313-1316 byatom
atomicfluct out substrate_rmsf_by_residue.dat :1313-1316 byres
distance d1_d_OG1_scissile_C :929@OG1 :1314@C out geometry_distances.dat
distance d2_d_Phe113N_scissile_O :113@N :1314@O out geometry_distances.dat
distance d3_d_Asn186ND2_scissile_O :186@ND2 :1314@O out geometry_distances.dat
angle a1_alpha_Burgi_Dunitz :929@OG1 :1314@C :1314@O out geometry_angles.dat
dihedral Herschlag_phi :1314@C5 :1315@N :1314@C :1314@O out geometry_angles.dat
rms ligand_rmsd :1313-1316&!@H= first out ligand_rmsd.dat
nativecontacts name protein_ligand_contacts :1-1312 :1313-1316 distance 4.0 out protein_ligand_contacts.dat
run
"""


def _minimization(restraint: str | None) -> str:
    title = "Restrained solute" if restraint else "Unrestrained"
    extra = ""
    if restraint:
        extra = f"  ntr=1, restraint_wt=10.0,\n  restraintmask='{restraint}',\n"
    return f"""{title} minimization
&cntrl
  imin=1, maxcyc=10000, ncyc=5000,
  ntb=1, cut=10.0,
{extra}/
"""


def _heating(seed: int, restraint: str) -> str:
    return f"""Restrained NVT heating 0 to 300 K
&cntrl
  imin=0, ntx=1, irest=0, nstlim=50000, dt=0.002,
  ntc=2, ntf=2, tempi=0.0, temp0=300.0,
  ntt=3, gamma_ln=1.0, ig={seed}, ntb=1, ntp=0, cut=10.0,
  ntpr=1000, ntwx=1000, ntwr=50000, ioutfm=1,
  ntr=1, restraint_wt=10.0, restraintmask='{restraint}',
/
&wt type='TEMP0', istep1=0, istep2=50000, value1=0.0, value2=300.0 /
&wt type='END' /
"""


def _equilibration(label: str, weight: float | None, restraint: str | None) -> str:
    extra = ""
    if weight is not None and restraint is not None:
        extra = f"  ntr=1, restraint_wt={weight:.1f}, restraintmask='{restraint}',\n"
    return f"""{label.capitalize()} NPT equilibration
&cntrl
  imin=0, ntx=5, irest=1, nstlim=500000, dt=0.002,
  ntc=2, ntf=2, temp0=300.0, ntt=3, gamma_ln=1.0,
  ntb=2, ntp=1, taup=2.0, cut=10.0,
  ntpr=5000, ntwx=5000, ntwr=500000, ioutfm=1,
{extra}/
"""


def _production(seed: int) -> str:
    return f"""Unrestrained 100 ns NPT production
&cntrl
  imin=0, ntx=5, irest=1, nstlim=50000000, dt=0.002,
  ntc=2, ntf=2, temp0=300.0, ntt=3, gamma_ln=1.0, ig={seed},
  ntb=2, ntp=1, taup=2.0, cut=10.0,
  ntpr=5000, ntwx=25000, ntwr=500000, ioutfm=1, ntxo=2,
/
"""
