"""Unrestrained Nyl12 E-S relaxation and Asp272 relay path utilities."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm


ATOMS = {
    "N": 13905,
    "H1": 13906,
    "H2": 13907,
    "OG1": 13916,
    "HG1": 13917,
    "ASP_OD1": 14559,
    "ASP_OD2": 14560,
    "C": 19903,
    "O": 19904,
    "AMIDE_N": 19905,
}
PHASE2_CORRECTION_REASON = (
    "Invalid phase II archived; restart required with retained relay coordinate."
)


def sha256(path: Path) -> str:
    """Return the SHA-256 identity of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    """Write finite, human-readable JSON metadata."""
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def distances(topology: Path, coordinates: Path) -> dict[str, float]:
    """Measure the reactive distances from one Amber coordinate set."""
    parm = AmberParm(str(topology), xyz=str(coordinates))
    xyz = parm.coordinates
    pairs = {
        "OG1_C": ("OG1", "C"),
        "OG1_HG1": ("OG1", "HG1"),
        "Asp272_OD1_HG1": ("ASP_OD1", "HG1"),
        "Asp272_OD2_HG1": ("ASP_OD2", "HG1"),
        "Thr_N_HG1": ("N", "HG1"),
        "Thr_N_H1": ("N", "H1"),
        "Thr_N_H2": ("N", "H2"),
        "C_O": ("C", "O"),
        "C_N": ("C", "AMIDE_N"),
        "OG1_Asp272_OD1": ("OG1", "ASP_OD1"),
    }
    return {
        name: float(np.linalg.norm(xyz[ATOMS[first] - 1] - xyz[ATOMS[second] - 1]))
        for name, (first, second) in pairs.items()
    }


def coordinates_valid(
    topology: Path, coordinates: Path, expected_box: list[float]
) -> bool:
    """Return whether all coordinates are finite and the fixed cell is unchanged."""
    parm = AmberParm(str(topology), xyz=str(coordinates))
    return bool(
        np.isfinite(parm.coordinates).all()
        and np.allclose(parm.box, expected_box, rtol=0, atol=1e-10)
    )


def es_input(base: str, max_cycles: int = 5000) -> str:
    """Create a fully unrestrained continuation input from the audited input."""
    if max_cycles < 1:
        raise ValueError("max_cycles must be positive")
    control = base.split("&dlfind", maxsplit=1)[0].rstrip()
    control = re.sub(r"\bntwx\s*=\s*\d+", "ntwx=25", control)
    text = (
        "Fully unrestrained noncovalent E-S continuation\n"
        + control[control.index("&cntrl") :]
        + "\n\n&dlfind optalg='LBFGS', crdrep='CART', active='@*', "
        + f"tol=1.d-4, tole=1.d-7, maxcycle={max_cycles}, "
        + "maxstep=0.5, verbosity=1, /\n"
    )
    validate_no_positional_restraints(text)
    return text


def driven_input(base: str, max_cycles: int = 2000) -> str:
    """Create a reaction-biased input with no positional restraints."""
    text = es_input(base, max_cycles=max_cycles)
    text, count = re.subn(r"\bnmropt\s*=\s*0", "nmropt=1", text, count=1)
    if count != 1:
        raise ValueError("Expected nmropt=0 in unrestrained base input")
    control, dlfind = text.split("&dlfind", maxsplit=1)
    text = (
        control.rstrip()
        + "\n&wt type='END' /\n"
        + "DISANG=restraints.RST\n"
        + "DUMPAVE=restraint_values.dat\n"
        + "&dlfind"
        + dlfind
    )
    validate_no_positional_restraints(text)
    return text


def validate_no_positional_restraints(text: str) -> None:
    """Fail if an input contains positional restraints or frozen atoms."""
    lowered = text.lower()
    forbidden = ("restraint_wt", "restraintmask", "ibelly", "bellymask")
    if any(token in lowered for token in forbidden):
        raise ValueError("Positional restraint or frozen-atom control found")
    if not re.search(r"\bntr\s*=\s*0\b", lowered):
        raise ValueError("Input must explicitly set ntr=0")


def schedule(initial: dict[str, float], windows: int = 21) -> list[dict[str, Any]]:
    """Return two bounded stages without extrapolating chemical bond targets."""
    if windows < 2:
        raise ValueError("At least two windows are required per stage")
    q1_start = initial["OG1_HG1"] - initial["Asp272_OD1_HG1"]
    rows: list[dict[str, Any]] = []
    for index in range(windows):
        fraction = index / (windows - 1)
        rows.append(
            {
                "label": f"stage1_{index:02d}",
                "stage": 1,
                "attack_target_A": initial["OG1_C"]
                + fraction * (1.55 - initial["OG1_C"]),
                "proton_difference_target_A": q1_start + fraction * (0.8 - q1_start),
            }
        )
    return rows


def stage2_schedule(
    initial_q2: float,
    relay_retention_target: float,
    og1_hg1_minimum: float,
    windows: int = 21,
) -> list[dict[str, Any]]:
    """Retain the first proton transfer while driving Asp272-H to Thr N."""
    if windows < 2:
        raise ValueError("At least two windows are required per stage")
    if not 0.3 <= relay_retention_target <= 2.0:
        raise ValueError("Relay-retention target must be between 0.3 and 2.0 A")
    if not 1.2 <= og1_hg1_minimum <= 2.0:
        raise ValueError("OG1-HG1 minimum must be between 1.2 and 2.0 A")
    return [
        {
            "label": f"stage2_{index:02d}",
            "stage": 2,
            "attack_target_A": 1.55,
            "proton_difference_target_A": initial_q2
            + index / (windows - 1) * (0.8 - initial_q2),
            "relay_retention_target_A": relay_retention_target,
            "og1_hg1_minimum_A": og1_hg1_minimum,
        }
        for index in range(windows)
    ]


def restraint_text(
    stage: int,
    attack: float,
    proton: float,
    relay_retention_target: float | None,
) -> str:
    """Build the declared harmonic reaction-coordinate restraints."""
    if stage not in {1, 2}:
        raise ValueError("stage must be 1 or 2")
    if stage == 1 and relay_retention_target is not None:
        raise ValueError("Stage I must not declare a relay-retention target")
    if stage == 2 and relay_retention_target is None:
        raise ValueError("Stage II requires a relay-retention target")
    proton_atoms = (
        f"{ATOMS['OG1']},{ATOMS['HG1']},{ATOMS['ASP_OD1']},{ATOMS['HG1']}"
        if stage == 1
        else f"{ATOMS['ASP_OD1']},{ATOMS['HG1']},{ATOMS['N']},{ATOMS['HG1']}"
    )
    blocks = []
    for atoms, weights, target in (
        (f"{ATOMS['OG1']},{ATOMS['C']},0,0,0,0,0,0", "1,0,0,0", attack),
        (proton_atoms + ",0,0,0,0", "1,-1,0,0", proton),
    ):
        blocks.append(
            f"&rst iat={atoms}, rstwt={weights}, r1={target - 2:.10f}, "
            f"r2={target:.10f}, r3={target:.10f}, r4={target + 2:.10f}, "
            "rk2=100, rk3=100, /\n"
        )
    if relay_retention_target is not None:
        atoms = f"{ATOMS['OG1']},{ATOMS['HG1']},{ATOMS['ASP_OD1']},{ATOMS['HG1']}"
        blocks.append(
            f"&rst iat={atoms},0,0,0,0, rstwt=1,-1,0,0, "
            f"r1={relay_retention_target - 2:.10f}, "
            f"r2={relay_retention_target:.10f}, "
            f"r3={relay_retention_target:.10f}, "
            f"r4={relay_retention_target + 2:.10f}, "
            "rk2=100, rk3=100, /\n"
        )
    return "".join(blocks)


def native_converged(output: Path) -> bool:
    """Return whether DL-FIND reports native convergence without a warning."""
    text = output.read_text()
    return (
        "Optimisation converged" in text
        and "WARNING: OPTIMISATION HAS NOT CONVERGED" not in text
    )
