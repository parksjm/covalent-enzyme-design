"""Match overlapping angular regions between relaxed scan directions."""

from __future__ import annotations

from collections.abc import Sequence


def match_profile_overlap(
    forward: Sequence[tuple[int, float]],
    reverse: Sequence[tuple[int, float]],
    forward_energies: Sequence[float],
    reverse_energies: Sequence[float],
    tolerance_degrees: float,
) -> list[dict[str, float | int]]:
    """Pair unique forward/reverse points within circular angle tolerance.

    Pair selection is deterministic: candidate pairs are sorted by circular
    angular distance and then by their input indices. Energies are compared
    only for matched angular overlap; unpaired path regions remain unassessed.

    Args:
        forward: Ordered ``(point_index, signed_angle_degrees)`` pairs.
        reverse: Ordered ``(point_index, signed_angle_degrees)`` pairs.
        forward_energies: Relative energies corresponding to forward points.
        reverse_energies: Relative energies corresponding to reverse points.
        tolerance_degrees: Maximum shortest circular separation in degrees.

    Returns:
        Matched points with signed angle and forward-minus-reverse energy.

    Raises:
        ValueError: If point/energy arrays do not align or tolerance is invalid.
    """
    if len(forward) != len(forward_energies) or len(reverse) != len(reverse_energies):
        raise ValueError("Scan point and energy arrays must have matching lengths")
    if tolerance_degrees < 0.0 or tolerance_degrees > 180.0:
        raise ValueError("Angular match tolerance must be between 0 and 180 degrees")
    candidates = [
        (_circular_distance(f_angle, r_angle), f_index, r_index)
        for f_index, (_, f_angle) in enumerate(forward)
        for r_index, (_, r_angle) in enumerate(reverse)
        if _circular_distance(f_angle, r_angle) <= tolerance_degrees
    ]
    used_forward: set[int] = set()
    used_reverse: set[int] = set()
    matches: list[dict[str, float | int]] = []
    for _, f_index, r_index in sorted(candidates):
        if f_index in used_forward or r_index in used_reverse:
            continue
        used_forward.add(f_index)
        used_reverse.add(r_index)
        matches.append(
            {
                "forward_point": forward[f_index][0],
                "reverse_point": reverse[r_index][0],
                "angle_degrees": forward[f_index][1],
                "hysteresis_kcal_mol": (
                    forward_energies[f_index] - reverse_energies[r_index]
                ),
            }
        )
    return sorted(matches, key=lambda row: int(row["forward_point"]))


def _circular_distance(first: float, second: float) -> float:
    """Return the shortest unsigned angular distance in degrees."""
    return abs((first - second + 180.0) % 360.0 - 180.0)
