"""Measure declared reaction distances without changing path coordinates."""

from __future__ import annotations

import math
from itertools import pairwise
from typing import Any

import numpy as np


def measure_distances(
    coordinates: np.ndarray, pairs: dict[str, dict[str, Any]]
) -> dict[str, float]:
    """Measure unwrapped Å distances for declared one-based atom pairs.

    Args:
        coordinates: Already imaged path coordinates, shape (atoms, 3).
        pairs: Named pairs with one-based ``atoms`` and descriptive groups.

    Returns:
        Distances in Å, without alignment or minimum-image transformation.
    """
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError("Invalid reaction coordinates")
    if not pairs:
        raise ValueError("Declare reaction distance pairs explicitly")
    result = {}
    for name, definition in pairs.items():
        a, b = definition["atoms"]
        if min(a, b) < 1 or max(a, b) > len(xyz) or a == b:
            raise ValueError(f"Invalid atom pair: {name}")
        result[name] = float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))
    return result


def distance_jumps(
    rows: list[dict[str, Any]], threshold: float
) -> list[dict[str, Any]]:
    """Flag adjacent-window changes for scientific review, not rejection.

    Args:
        rows: Ordered endpoint records with window and distances_A keys.
        threshold: Strictly positive absolute distance-change threshold in Å.

    Returns:
        Changes exceeding the threshold, retaining stage boundaries.
    """
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("Jump threshold must be finite and positive")
    flags = []
    for previous, current in pairwise(rows):
        if previous["distances_A"].keys() != current["distances_A"].keys():
            raise ValueError("Distance definitions changed between windows")
        for name, value in current["distances_A"].items():
            delta = value - previous["distances_A"][name]
            if not math.isfinite(delta):
                raise ValueError("Nonfinite reaction distance change")
            if abs(delta) > threshold:
                flags.append(
                    {
                        "from": previous["window"],
                        "to": current["window"],
                        "distance": name,
                        "delta_A": delta,
                        "stage_boundary": previous["stage"] != current["stage"],
                    }
                )
    return flags
