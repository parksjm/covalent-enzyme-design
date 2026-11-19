"""Mapped direct OG1-to-substrate-N acyl-enzyme diagnostic utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from parmed.amber import AmberParm

from nyl12_covalent.coupled_ti1_drive import coupled_geometry

SUBSTRATE_N = 19929
HG1 = 13917


def direct_geometry(topology: Path, coordinates: Path) -> dict[str, float]:
    """Measure reaction pairs and retained QM hydrogen-parent identities."""
    values = coupled_geometry(topology, coordinates)
    system = AmberParm(str(topology), xyz=str(coordinates))
    xyz = system.coordinates

    def distance(a: int, b: int) -> float:
        return float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))

    nitrogen = system.atoms[SUBSTRATE_N - 1]
    hydrogens = [a for a in nitrogen.bond_partners if a.atomic_number == 1]
    if len(hydrogens) != 1 or hydrogens[0].idx + 1 != 19937:
        raise ValueError("Expected substrate N19929 original HN19937")
    values["substrate_N_HG1_A"] = distance(SUBSTRATE_N, HG1)
    values["substrate_N_HN_A"] = distance(SUBSTRATE_N, 19937)
    pairs = []
    qm = [
        *range(1621, 1656),
        *range(2727, 2735),
        *range(13905, 13932),
        *range(14555, 14561),
        *range(14574, 14580),
        *range(19914, 19918),
        *range(19925, 19932),
        *range(19937, 19942),
    ]
    for index in qm:
        atom = system.atoms[index - 1]
        if atom.atomic_number != 1 or index == HG1:
            continue
        if len(atom.bond_partners) != 1:
            raise ValueError(f"Expected one hydrogen parent for atom {index}")
        pairs.append(distance(index, atom.bond_partners[0].idx + 1))
    values["retained_QM_H_parent_min_A"] = min(pairs)
    values["retained_QM_H_parent_max_A"] = max(pairs)
    values["proton_difference_A"] = values["OG1_HG1_A"] - values["substrate_N_HG1_A"]
    return values


def composite(values: dict[str, float]) -> float:
    """Return proton transfer plus formation and cleavage coordinate."""
    return (
        values["proton_difference_A"]
        - values["OG1_carbonyl_C_A"]
        + values["carbonyl_C_N_A"]
    )


def chemistry_retained(values: dict[str, float]) -> bool:
    """Allow mapped asynchronous chemistry while retaining other QM protons."""
    return (
        0.75 <= values["retained_QM_H_parent_min_A"]
        and values["retained_QM_H_parent_max_A"] <= 1.30
        and all(0.85 <= values[f"N_{h}_A"] <= 1.25 for h in ("H1", "H2"))
        and 0.85 <= values["substrate_N_HN_A"] <= 1.25
        and values["N_HG1_A"] >= 1.30
        and all(values[f"OG1_{h}_A"] >= 1.30 for h in ("H1", "H2"))
        and 0.75 <= min(values["substrate_N_HG1_A"], values["OG1_HG1_A"]) <= 1.70
        and 1.25 <= values["OG1_carbonyl_C_A"] <= 4.00
        and 1.10 <= values["carbonyl_C_O_A"] <= 1.55
        and 1.20 <= values["carbonyl_C_N_A"] <= 4.50
    )


def candidate_geometry(values: dict[str, float]) -> bool:
    """Identify an ester/NH2 candidate for complete reaction-bias release."""
    return (
        chemistry_retained(values)
        and values["OG1_carbonyl_C_A"] <= 1.60
        and values["carbonyl_C_N_A"] >= 2.40
        and values["carbonyl_C_O_A"] <= 1.30
        and 0.85 <= values["substrate_N_HG1_A"] <= 1.25
        and values["OG1_HG1_A"] >= 1.50
    )


def direct_restraint(
    target: float,
    force_constant: float = 100.0,
    *,
    substrate_proton_guard: bool = False,
) -> str:
    """Restrain only the four-pair asynchronous direct-acyl coordinate."""
    if force_constant <= 0:
        raise ValueError("force_constant must be positive")
    text = (
        "&rst iat=13916,13917,13917,19929,13916,19916,19916,19929, "
        "rstwt=1,-1,-1,1, "
        f"r1={target - 2:.10f}, r2={target:.10f}, r3={target:.10f}, "
        f"r4={target + 2:.10f}, rk2={force_constant:g}, "
        f"rk3={force_constant:g}, /\n"
    )
    if substrate_proton_guard:
        text += (
            "&rst iat=19929,13917, r1=0.8, r2=0.95, r3=1.25, r4=1.5, "
            "rk2=50, rk3=50, /\n"
        )
    return text
