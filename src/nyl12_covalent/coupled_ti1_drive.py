"""Coupled direct proton-transfer and nucleophilic-attack path utilities."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm

from nyl12_covalent.direct_proton_transfer import (
    HG1_ATOM,
    N_ATOM,
    OG1_ATOM,
    proton_distances,
)

CARBONYL_C_ATOM = 19916
CARBONYL_O_ATOM = 19917
H1_ATOM = 13906
H2_ATOM = 13907


def coupled_geometry(topology: Path, coordinates: Path) -> dict[str, float]:
    """Measure the coupled proton-transfer and nucleophilic-attack coordinates."""
    values = proton_distances(topology, coordinates)
    system = AmberParm(str(topology), xyz=str(coordinates))
    xyz = system.coordinates
    values["OG1_carbonyl_C_A"] = float(
        np.linalg.norm(xyz[OG1_ATOM - 1] - xyz[CARBONYL_C_ATOM - 1])
    )
    values["carbonyl_C_O_A"] = float(
        np.linalg.norm(xyz[CARBONYL_C_ATOM - 1] - xyz[CARBONYL_O_ATOM - 1])
    )
    carbon = system.atoms[CARBONYL_C_ATOM - 1]
    nitrogens = [atom for atom in carbon.bond_partners if atom.atomic_number == 7]
    if len(nitrogens) != 1:
        raise ValueError("Expected one scissile amide nitrogen bonded to carbonyl C")
    values["carbonyl_C_N_A"] = float(
        np.linalg.norm(xyz[carbon.idx] - xyz[nitrogens[0].idx])
    )
    for label, atom in (("H1", H1_ATOM), ("H2", H2_ATOM)):
        values[f"N_{label}_A"] = float(np.linalg.norm(xyz[N_ATOM - 1] - xyz[atom - 1]))
        values[f"OG1_{label}_A"] = float(
            np.linalg.norm(xyz[OG1_ATOM - 1] - xyz[atom - 1])
        )
    return values


def ti1_like_geometry(values: dict[str, float]) -> bool:
    """Return whether driven coordinates have reached a TI1-like endpoint."""
    return (
        values["OG1_carbonyl_C_A"] <= 1.80
        and values["N_HG1_A"] <= 1.25
        and values["OG1_HG1_A"] >= 1.30
        and values["carbonyl_C_O_A"] >= 1.28
    )


def coupled_schedule(
    initial_proton: float,
    initial_attack: float,
    final_proton: float = 0.75,
    final_attack: float = 1.55,
    maximum_increment: float = 0.10,
) -> list[dict[str, Any]]:
    """Build a synchronous path with small increments in both coordinates."""
    if final_proton <= initial_proton or final_attack >= initial_attack:
        raise ValueError("Coupled path endpoints must advance toward TI1")
    if maximum_increment <= 0.0:
        raise ValueError("maximum_increment must be positive")
    intervals = max(
        math.ceil((final_proton - initial_proton) / maximum_increment),
        math.ceil((initial_attack - final_attack) / maximum_increment),
    )
    return [
        {
            "label": f"coupled_{index:02d}",
            "proton_target_A": initial_proton
            + index / intervals * (final_proton - initial_proton),
            "attack_target_A": initial_attack
            + index / intervals * (final_attack - initial_attack),
        }
        for index in range(1, intervals + 1)
    ]


def coupled_restraints(
    proton_target: float,
    attack_target: float,
    force_constant: float = 100.0,
) -> str:
    """Return harmonic restraints for the two declared reaction coordinates."""
    if force_constant <= 0.0:
        raise ValueError("force_constant must be positive")
    proton = (
        f"&rst iat={OG1_ATOM},{HG1_ATOM},{N_ATOM},{HG1_ATOM},0,0,0,0, "
        f"rstwt=1,-1,0,0, r1={proton_target - 2.0:.10f}, "
        f"r2={proton_target:.10f}, r3={proton_target:.10f}, "
        f"r4={proton_target + 2.0:.10f}, rk2={force_constant:g}, "
        f"rk3={force_constant:g}, /\n"
    )
    attack = (
        f"&rst iat={OG1_ATOM},{CARBONYL_C_ATOM},0,0,0,0,0,0, "
        f"rstwt=1,0,0,0, r1={attack_target - 2.0:.10f}, "
        f"r2={attack_target:.10f}, r3={attack_target:.10f}, "
        f"r4={attack_target + 2.0:.10f}, rk2={force_constant:g}, "
        f"rk3={force_constant:g}, /\n"
    )
    return proton + attack


def asynchronous_schedule(
    initial_proton: float,
    initial_attack: float,
    final_composite: float = -0.80,
    maximum_increment: float = 0.15,
) -> list[dict[str, Any]]:
    """Drive one composite coordinate while component distances remain free."""
    initial = initial_proton - initial_attack
    if final_composite <= initial or maximum_increment <= 0.0:
        raise ValueError("Asynchronous endpoint must advance with a positive step")
    intervals = math.ceil((final_composite - initial) / maximum_increment)
    return [
        {
            "label": f"asynchronous_{index:02d}",
            "composite_target_A": initial
            + index / intervals * (final_composite - initial),
        }
        for index in range(1, intervals + 1)
    ]


def asynchronous_extension_schedule(
    last_target: float,
    target_increment: float,
    first_index: int,
    maximum_windows: int = 8,
) -> list[dict[str, Any]]:
    """Extend an asynchronous path with the established target spacing.

    Args:
        last_target: Final composite-coordinate target from the original path.
        target_increment: Positive spacing used by the original path.
        first_index: Sequential label number for the first extension window.
        maximum_windows: Safety limit; execution may stop earlier at TI1 geometry.

    Returns:
        Ordered extension targets with continuous asynchronous labels.
    """
    if target_increment <= 0.0:
        raise ValueError("target_increment must be positive")
    if first_index < 1 or maximum_windows < 1:
        raise ValueError("Positive label and window limits are required")
    return [
        {
            "label": f"asynchronous_{first_index + offset:02d}",
            "composite_target_A": last_target + (offset + 1) * target_increment,
        }
        for offset in range(maximum_windows)
    ]


def asynchronous_restraint(
    composite_target: float, force_constant: float = 100.0
) -> str:
    """Restrain only q_PT minus the OG1–carbonyl attack distance."""
    if force_constant <= 0.0:
        raise ValueError("force_constant must be positive")
    return (
        f"&rst iat={OG1_ATOM},{HG1_ATOM},{N_ATOM},{HG1_ATOM},"
        f"{OG1_ATOM},{CARBONYL_C_ATOM},0,0, rstwt=1,-1,-1,0, "
        f"r1={composite_target - 2.0:.10f}, r2={composite_target:.10f}, "
        f"r3={composite_target:.10f}, r4={composite_target + 2.0:.10f}, "
        f"rk2={force_constant:g}, rk3={force_constant:g}, /\n"
    )


def attack_only_schedule(
    initial_attack: float,
    final_attack: float = 1.50,
    maximum_increment: float = 0.075,
) -> list[dict[str, Any]]:
    """Build a bounded attack schedule from a proton-transferred geometry.

    Args:
        initial_attack: Initial OG1--carbonyl-C distance in angstrom.
        final_attack: Final target distance in angstrom.
        maximum_increment: Largest decrease between adjacent targets.

    Returns:
        Ordered target labels and distances.
    """
    if final_attack >= initial_attack or maximum_increment <= 0.0:
        raise ValueError("Attack endpoint must be shorter with a positive step")
    intervals = math.ceil((initial_attack - final_attack) / maximum_increment)
    return [
        {
            "label": f"attack_{index:02d}",
            "attack_target_A": initial_attack
            + index / intervals * (final_attack - initial_attack),
        }
        for index in range(1, intervals + 1)
    ]


def attack_with_proton_identity_restraints(
    attack_target: float,
    attack_force_constant: float = 100.0,
    proton_force_constant: float = 50.0,
) -> str:
    """Drive OG1 attack while retaining all three terminal N--H identities.

    The N--H restraints are flat-bottom gates, not fixed bond lengths. All
    three protons move freely between 0.90 and 1.20 A from the nitrogen.

    Args:
        attack_target: Harmonic OG1--carbonyl-C target in angstrom.
        attack_force_constant: Attack restraint force constant.
        proton_force_constant: N--H flat-bottom wall force constant.

    Returns:
        Amber NMR restraint records.
    """
    if attack_force_constant <= 0.0 or proton_force_constant <= 0.0:
        raise ValueError("Restraint force constants must be positive")
    records = [
        f"&rst iat={OG1_ATOM},{CARBONYL_C_ATOM},0,0,0,0,0,0, "
        f"rstwt=1,0,0,0, r1={attack_target - 2.0:.10f}, "
        f"r2={attack_target:.10f}, r3={attack_target:.10f}, "
        f"r4={attack_target + 2.0:.10f}, rk2={attack_force_constant:g}, "
        f"rk3={attack_force_constant:g}, /"
    ]
    for hydrogen in (HG1_ATOM, H1_ATOM, H2_ATOM):
        records.append(
            f"&rst iat={N_ATOM},{hydrogen},0,0,0,0,0,0, "
            "rstwt=1,0,0,0, r1=0.8000000000, r2=0.9000000000, "
            "r3=1.2000000000, r4=1.3500000000, "
            f"rk2={proton_force_constant:g}, rk3={proton_force_constant:g}, /"
        )
    return "\n".join(records) + "\n"
