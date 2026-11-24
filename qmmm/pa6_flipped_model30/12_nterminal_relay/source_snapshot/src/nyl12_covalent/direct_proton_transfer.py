"""Direct Thr233 O-gamma-to-N-terminal proton-transfer utilities."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm

N_ATOM = 13905
OG1_ATOM = 13916
HG1_ATOM = 13917


def proton_distances(topology: Path, coordinates: Path) -> dict[str, float]:
    """Measure the direct proton-transfer distances and difference coordinate.

    Args:
        topology: Amber topology for the solvated system.
        coordinates: Amber restart containing one coordinate frame.

    Returns:
        Distances in angstrom and ``d(OG1,HG1) - d(N,HG1)``.
    """
    system = AmberParm(str(topology), xyz=str(coordinates))
    xyz = system.coordinates
    donor = float(np.linalg.norm(xyz[OG1_ATOM - 1] - xyz[HG1_ATOM - 1]))
    acceptor = float(np.linalg.norm(xyz[N_ATOM - 1] - xyz[HG1_ATOM - 1]))
    return {
        "OG1_HG1_A": donor,
        "N_HG1_A": acceptor,
        "proton_difference_A": donor - acceptor,
    }


def schedule(
    initial_difference: float,
    final_difference: float = 0.75,
    maximum_increment: float = 0.10,
) -> list[dict[str, Any]]:
    """Build a bounded direct proton-transfer schedule.

    Args:
        initial_difference: Initial difference coordinate in angstrom.
        final_difference: Product-side difference coordinate in angstrom.
        maximum_increment: Largest allowed target increment in angstrom.

    Returns:
        Sequential window labels and difference-coordinate targets. The first
        target advances from the starting structure rather than duplicating it.
    """
    if final_difference <= initial_difference:
        raise ValueError("Final proton coordinate must exceed the initial value")
    if maximum_increment <= 0.0:
        raise ValueError("maximum_increment must be positive")
    intervals = math.ceil((final_difference - initial_difference) / maximum_increment)
    return [
        {
            "label": f"pt_{index:02d}",
            "target_A": initial_difference
            + index / intervals * (final_difference - initial_difference),
        }
        for index in range(1, intervals + 1)
    ]


def restraint_text(target: float, force_constant: float = 100.0) -> str:
    """Return one generalized distance-difference restraint.

    Args:
        target: Target for ``d(OG1,HG1) - d(N,HG1)`` in angstrom.
        force_constant: Flat-bottom restraint force constant.

    Returns:
        Amber NMR-style restraint block.
    """
    if force_constant <= 0.0:
        raise ValueError("force_constant must be positive")
    return (
        f"&rst iat={OG1_ATOM},{HG1_ATOM},{N_ATOM},{HG1_ATOM},0,0,0,0, "
        f"rstwt=1,-1,0,0, r1={target - 2.0:.10f}, r2={target:.10f}, "
        f"r3={target:.10f}, r4={target + 2.0:.10f}, "
        f"rk2={force_constant:g}, rk3={force_constant:g}, /\n"
    )


def driven_input(base: str, max_cycles: int = 2_000) -> str:
    """Enable only the generalized proton-transfer restraint in an audited input.

    Args:
        base: Accepted 12-angstrom-shell DFTB3/MM input.
        max_cycles: Maximum DL-FIND cycles per window.

    Returns:
        Driven DFTB3/MM input retaining the original QM and active-region masks.
    """
    _validate_base(base)
    text = _set_cycle_limits(base, max_cycles)
    text, count = re.subn(r"\bnmropt\s*=\s*0", "nmropt=1", text, count=1)
    if count != 1:
        raise ValueError("Expected exactly one nmropt=0 control")
    control, dlfind = text.split("&dlfind", maxsplit=1)
    return (
        control.rstrip()
        + "\n&wt type='END' /\nDISANG=restraints.RST\n"
        + "DUMPAVE=restraint_values.dat\n&dlfind"
        + dlfind
    )


def release_input(base: str, max_cycles: int = 10_000) -> str:
    """Return the accepted 12-angstrom-shell input with no reaction bias.

    Args:
        base: Accepted 12-angstrom-shell DFTB3/MM input.
        max_cycles: Maximum cycles for the endpoint relaxation.

    Returns:
        Unbiased input retaining the original QM and active-region masks.
    """
    _validate_base(base)
    return _set_cycle_limits(base, max_cycles)


def _set_cycle_limits(text: str, max_cycles: int) -> str:
    if max_cycles < 1:
        raise ValueError("max_cycles must be positive")
    result, cntrl_count = re.subn(
        r"\bmaxcyc\s*=\s*\d+", f"maxcyc={max_cycles}", text, count=1
    )
    result, dlfind_count = re.subn(
        r"\bmaxcycle\s*=\s*\d+", f"maxcycle={max_cycles}", result, count=1
    )
    if (cntrl_count, dlfind_count) != (1, 1):
        raise ValueError("Expected one maxcyc and one maxcycle control")
    return result


def _validate_base(text: str) -> None:
    lowered = text.lower()
    required = ("nmropt=0", "ntr=0", "qm_theory='dftb3'", "active='")
    if any(token not in lowered.replace(" ", "") for token in required[:3]):
        raise ValueError("Base input is not the expected unbiased DFTB3/MM input")
    if required[3] not in lowered:
        raise ValueError("Base input lacks a declared DL-FIND active region")
    forbidden = ("restraint_wt", "restraintmask", "ibelly", "bellymask", "disang")
    if any(token in lowered for token in forbidden):
        raise ValueError("Base input contains an undeclared restraint or frozen mask")
