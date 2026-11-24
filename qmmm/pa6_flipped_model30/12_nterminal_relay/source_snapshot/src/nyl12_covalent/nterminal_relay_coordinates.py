"""Coordinates and conservative endpoint gates for N-terminal relay exploration.

These gates identify geometric candidates only; they do not establish a minimum,
transition state, complete chemical graph, or valid proton inventory elsewhere.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from nyl12_covalent.relay_release import ROLES, geometry as relay_geometry

_BRANCH_ACCEPTOR = {"reverse_activation": "OG", "forward_acylation": "NL"}
_PROTONS = (ROLES["H1"], ROLES["H2"], ROLES["HG"])


def _acceptor(branch: str, proton: int) -> int:
    if branch not in _BRANCH_ACCEPTOR:
        raise ValueError(f"Unknown relay branch: {branch}")
    if proton not in _PROTONS:
        raise ValueError("Proton must be one of the three original Nalpha hydrogens")
    return ROLES[_BRANCH_ACCEPTOR[branch]]


def coordinate(xyz: np.ndarray, branch: str, proton: int = 13917) -> float:
    """Calculate the composite transfer and bond distance in angstrom.

    Args:
        xyz: Full original-topology coordinates in angstrom.
        branch: Either reverse_activation or forward_acylation.
        proton: Explicit one-based original Nalpha hydrogen index.

    Returns:
        Nalpha–H minus acceptor–H plus carbon–acceptor distance.
    """
    acceptor = _acceptor(branch, proton)
    if xyz.shape != (104922, 3):
        raise ValueError("Full original-topology coordinates required")
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite coordinates")
    pos = lambda i: xyz[i - 1]  # noqa: E731
    return float(
        np.linalg.norm(pos(ROLES["N"]) - pos(proton))
        - np.linalg.norm(pos(acceptor) - pos(proton))
        + np.linalg.norm(pos(ROLES["C"]) - pos(acceptor))
    )


def amber_restraint(
    target: float, branch: str, proton: int = 13917, force: float = 100
) -> str:
    """Build an Amber generalized-distance restraint for one relay window.

    Args:
        target: Composite coordinate target in angstrom.
        branch: Either reverse_activation or forward_acylation.
        proton: Explicit one-based original Nalpha hydrogen index.
        force: Positive Amber restraint force constant.

    Returns:
        One complete Amber DISANG restraint record.
    """
    acceptor = _acceptor(branch, proton)
    if not math.isfinite(target) or not math.isfinite(force) or force <= 0:
        raise ValueError("Target must be finite and force positive and finite")
    return (
        f"&rst iat={ROLES['N']},{proton},{acceptor},{proton},"
        f"{ROLES['C']},{acceptor},0,0, rstwt=1,-1,1,0, "
        f"r1={target - 2:.10f}, r2={target:.10f}, r3={target:.10f}, "
        f"r4={target + 2:.10f}, rk2={force:.10f}, rk3={force:.10f}, /\n"
    )


def schedule(
    initial: float, final: float = 4.2, maxincrement: float = 0.15
) -> list[dict[str, Any]]:
    """Generate ascending named windows, excluding the initial coordinate.

    Args:
        initial: Starting composite coordinate.
        final: Last composite coordinate, strictly greater than initial.
        maxincrement: Positive maximum spacing between consecutive targets.

    Returns:
        Ordered records containing window names and coordinate targets.
    """
    if not all(math.isfinite(v) for v in (initial, final, maxincrement)):
        raise ValueError("Schedule values must be finite")
    if final <= initial or maxincrement <= 0:
        raise ValueError("Require final > initial and positive maxincrement")
    count = math.ceil((final - initial) / maxincrement)
    return [
        {"window": f"window_{i:03d}", "target": initial + (final - initial) * i / count}
        for i in range(1, count + 1)
    ]


def geometry(xyz: np.ndarray) -> dict[str, Any]:
    """Measure existing relay geometry and conservative proton ownership.

    Args:
        xyz: Full original-topology coordinates in angstrom.

    Returns:
        Existing relay measurements with Nalpha–leaving-N distance and
        1.30 Å ownership lists for each of the four tracked hydrogens.
    """
    if xyz.shape != (104922, 3):
        raise ValueError("Full original-topology coordinates required")
    result = relay_geometry(xyz)
    result["distances_A"]["N_NL"] = float(
        np.linalg.norm(xyz[ROLES["N"] - 1] - xyz[ROLES["NL"] - 1])
    )
    distances = result["distances_A"]
    result["ownership_cutoff_A"] = 1.30
    result["ownership"] = {
        h: [a for a in ("N", "OG", "NL", "O") if distances[f"{h}_{a}"] <= 1.30]
        for h in ("H1", "H2", "HG", "HN")
    }
    return result


def candidate(geometry: dict[str, Any], branch: str) -> bool:
    """Classify a geometric endpoint without identifying which proton transfers.

    Args:
        geometry: Result from :func:`geometry`.
        branch: Either reverse_activation or forward_acylation.

    Returns:
        True if distances and unambiguous tracked proton ownership meet gates.
    """
    if branch not in _BRANCH_ACCEPTOR:
        raise ValueError(f"Unknown relay branch: {branch}")
    d, owners = geometry["distances_A"], geometry["ownership"]
    if any(len(owners[h]) != 1 for h in ("H1", "H2", "HG", "HN")):
        return False
    if owners["HN"] != ["NL"]:
        return False
    original = [owners[h][0] for h in ("H1", "H2", "HG")]
    acceptor = _BRANCH_ACCEPTOR[branch]
    if original.count("N") != 2 or original.count(acceptor) != 1:
        return False
    if not 1.1 <= d["C_O"] <= 1.32:
        return False
    if branch == "reverse_activation":
        return d["C_OG"] >= 2.5 and 1.2 <= d["C_NL"] <= 1.65
    return 1.25 <= d["C_OG"] <= 1.65 and d["C_NL"] >= 2.2
