"""Audited decomposition of the optional direct-acyl proton guard energy."""

from __future__ import annotations

import re


def proton_guard_bias(restraints: str, distance: float) -> float | None:
    """Evaluate an explicitly present, mapped substrate-N/HG1 flat-bottom bias.

    Args:
        restraints: Amber restraint namelists for one window.
        distance: Measured substrate N19929--HG113917 distance in angstroms.

    Returns:
        The guard energy in kcal/mol, or None when no guard is declared.

    Raises:
        ValueError: If the declared guard differs from the audited protocol.
    """
    blocks = re.findall(r"&rst\b(.*?)/", restraints, flags=re.S | re.I)
    guards = [
        block for block in blocks if re.search(r"iat\s*=\s*19929\s*,\s*13917", block)
    ]
    if not guards:
        return None
    if len(guards) != 1 or len(blocks) != 2:
        raise ValueError(
            "Expected exactly composite and substrate proton-guard restraints"
        )
    parameters = {
        key.lower(): float(value)
        for key, value in re.findall(r"\b(r[1-4]|rk[23])\s*=\s*([-+\d.eE]+)", guards[0])
    }
    if parameters != dict(r1=0.8, r2=0.95, r3=1.25, r4=1.5, rk2=50, rk3=50):
        raise ValueError("Substrate proton guard differs from audited protocol")
    if distance < 0.8:
        return 50 * (0.8 - 0.95) ** 2 + 100 * (0.8 - 0.95) * (distance - 0.8)
    if distance < 0.95:
        return 50 * (distance - 0.95) ** 2
    if distance <= 1.25:
        return 0.0
    if distance <= 1.5:
        return 50 * (distance - 1.25) ** 2
    return 50 * (1.5 - 1.25) ** 2 + 100 * (1.5 - 1.25) * (distance - 1.5)
