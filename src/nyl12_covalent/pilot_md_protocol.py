"""Matched 300 K Nyl12 MD pilot inputs for protonation comparisons."""

from __future__ import annotations

from typing import Any

STAGES = (
    ("01_heat_nvt", 50_000, 5.0),
    ("02_density_npt", 50_000, 2.0),
    ("03_weak_npt", 50_000, 0.5),
    ("04_equil_npt", 100_000, None),
    ("05_production_10ns", 5_000_000, None),
)
HEAVY_SOLUTE_MASK = "!:WAT,Na+,Cl- & !@H="


def protocol_manifest(seed: int) -> dict[str, Any]:
    """Return the common protocol and replica-specific random seeds.

    Args:
        seed: Independent positive heating seed.

    Returns:
        Explicit common protocol and all stage specifications.
    """
    if seed <= 0:
        raise ValueError("Random seed must be positive")
    return {
        "protocol": "model30_matched_10ns_pilot_v1",
        "dt_ps": 0.002,
        "temperature_K": 300,
        "thermostat": "Langevin ntt=3 gamma_ln=1.0 ps^-1",
        "barostat": "Monte Carlo barostat=2 ntp=1 pres0=1.0 taup=2.0 ps",
        "constraints": "SHAKE hydrogen bonds ntc=2 ntf=2",
        "electrostatics": "Periodic PME cut=10.0 A",
        "heavy_restraint_mask": HEAVY_SOLUTE_MASK,
        "restraint_reference": "Original minimized coordinates for every restrained stage",
        "production_restraints": "ntr=0 nmropt=0",
        "trajectory_interval_ps": {"equilibration": 2, "production": 10},
        "stages": [
            {
                "name": name,
                "steps": steps,
                "duration_ps": steps * 0.002,
                "heavy_restraint_kcal_mol_A2": weight,
                "seed": seed + 10_000 * index,
            }
            for index, (name, steps, weight) in enumerate(STAGES)
        ],
        "analysis": "Use projected_herschlag_phi, not old four-atom cpptraj dihedral",
    }


def mdin_files(seed: int) -> dict[str, str]:
    """Return100ps heating,400ps equilibration and10ns production inputs.

    Args:
        seed: Independent positive heating seed.

    Returns:
        Complete Amber input text keyed by stage.
    """
    inputs = {}
    for index, stage in enumerate(protocol_manifest(seed)["stages"]):
        heating = index == 0
        lines = [
            stage["name"],
            "&cntrl",
            " imin=0,",
            f" ntx={1 if heating else 5}, irest={0 if heating else 1},",
            f" nstlim={stage['steps']}, dt=0.002, ntc=2, ntf=2,",
            f" ntt=3, gamma_ln=1.0, temp0=300.0, ig={stage['seed']},",
            f" ntb={1 if heating else 2}, ntp={0 if heating else 1}, cut=10.0,",
            f" ntpr=1000, ntwx={5000 if index == 4 else 1000}, ntwr=50000,",
            " ioutfm=1, ntxo=2, iwrap=0,",
        ]
        lines.append(
            " tempi=0.0, nmropt=1,"
            if heating
            else " barostat=2, pres0=1.0, taup=2.0, nmropt=0,"
        )
        weight = stage["heavy_restraint_kcal_mol_A2"]
        lines.append(
            " ntr=0,"
            if weight is None
            else f" ntr=1, restraint_wt={weight}, restraintmask='{HEAVY_SOLUTE_MASK}',"
        )
        lines.append("/")
        if heating:
            lines.extend(
                (
                    "&wt type='TEMP0', istep1=0, istep2=50000, value1=0.0, value2=300.0 /",
                    "&wt type='END' /",
                )
            )
        inputs[stage["name"]] = "\n".join(lines) + "\n"
    return inputs
