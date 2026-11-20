"""Native Amber two-coordinate umbrellas and exact cross-window bias energies.

Lengths use Angstrom and energies kcal/mol. Force constants are Amber rk values:
the central potential is rk*(q-center)**2, without a factor of one half.
Outside center +/- outer_extent the potential has Amber's tangent linear tails.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping


@dataclass(frozen=True)
class HarmonicWindow:
    """Two independent native Amber umbrellas with zero-width flat bottoms."""

    q_acyl_center: float
    q_h_center: float
    k_acyl: float
    k_h: float
    outer_extent: float = 10.0

    def __post_init__(self) -> None:
        values = (
            self.q_acyl_center,
            self.q_h_center,
            self.k_acyl,
            self.k_h,
            self.outer_extent,
        )
        if not all(math.isfinite(v) for v in values):
            raise ValueError("Window values must be finite")
        if min(self.k_acyl, self.k_h, self.outer_extent) <= 0:
            raise ValueError("Amber rk constants and outer extent must be positive")


def bias_energy(cv_values: tuple[float, float], window: HarmonicWindow) -> float:
    """Return exact two-restraint native energy, including linear outer tails.

    Args:
        cv_values: (qAcyl, qH), both in Angstrom.
        window: Amber rk constants and centers, including outer-tail boundaries.

    Returns:
        Total bias energy in kcal/mol; directly usable for all-window MBAR.
    """
    if len(cv_values) != 2 or not all(math.isfinite(q) for q in cv_values):
        raise ValueError("Need two finite coordinates in Angstrom")
    total = 0.0
    for q, center, k in zip(
        cv_values,
        (window.q_acyl_center, window.q_h_center),
        (window.k_acyl, window.k_h),
        strict=True,
    ):
        d = abs(q - center)
        width = window.outer_extent
        total += k * (d * d if d <= width else width * (2 * d - width))
    return total


def amber_restraints(window: HarmonicWindow, roles: Mapping[str, int]) -> str:
    """Build qAcyl=C-Nleave minus C-OG1 and qH=OG1-H minus Nleave-H.

    Args:
        window: Explicit centers and Amber rk force constants.
        roles: One-based carbon, leaving_n, thr_og, proton atom IDs.

    Returns:
        Two independent generalized-distance &rst records.
    """
    keys = ("carbon", "leaving_n", "thr_og", "proton")
    atoms = [roles[k] for k in keys]
    if any(type(a) is not int or a <= 0 for a in atoms) or len(set(atoms)) != 4:
        raise ValueError("Four distinct positive one-based atom IDs required")
    c, n, o, h = atoms
    records = []
    for pairs, center, k in (
        ((c, n, c, o), window.q_acyl_center, window.k_acyl),
        ((o, h, n, h), window.q_h_center, window.k_h),
    ):
        ids = ",".join(str(a) for a in (*pairs, 0, 0, 0, 0))
        records.append(
            f"&rst iat={ids}, rstwt=1,-1,0,0, ialtd=0, "
            f"r1={center - window.outer_extent:.12g}, r2={center:.12g}, "
            f"r3={center:.12g}, r4={center + window.outer_extent:.12g}, "
            f"rk2={k:.12g}, rk3={k:.12g}, /"
        )
    return "\n".join(records) + "\n"
