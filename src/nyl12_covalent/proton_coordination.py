"""Smooth proton occupancy coordinates with analytic Cartesian derivatives.

These coordinates measure geometric coordination, not formal charge or pKa.
Hydrogen inputs must represent real atoms from the declared QM region.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class CoordinationTransfer:
    """Two smooth occupancies and their difference/gradient in inverse Å."""

    donor: float
    acceptor: float
    value: float
    gradient: FloatArray


def smooth_coordination(
    coordinates: FloatArray,
    center: int,
    hydrogens: Sequence[int],
    radius: float = 1.25,
) -> tuple[float, FloatArray]:
    """Evaluate sum_H 1/(1+(r_CH/radius)^6) and its exact derivative.

    Args:
        coordinates: Cartesian array with shape (N, 3), in Å.
        center: Zero-based heavy-atom index.
        hydrogens: Distinct zero-based real hydrogen indices, excluding center.
        radius: Positive switching length, in Å; counts are not integers.

    Returns:
        Smooth count and Cartesian derivative array with shape (N, 3).

    Raises:
        ValueError: If coordinates, indices, or switching radius are invalid.
    """
    xyz = np.asarray(coordinates, dtype=float)
    raw_h = np.asarray(hydrogens)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError("Coordinates must be a finite (N,3) array")
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Switching radius must be finite and positive")
    if not isinstance(center, (int, np.integer)) or not 0 <= center < len(xyz):
        raise ValueError("Center index is invalid")
    if raw_h.ndim != 1 or not len(raw_h) or not np.issubdtype(raw_h.dtype, np.integer):
        raise ValueError("Hydrogen indices must be a nonempty integer vector")
    h = raw_h.astype(int)
    if len(set(h)) != len(h) or np.any(h < 0) or np.any(h >= len(xyz)) or center in h:
        raise ValueError("Hydrogen indices must be distinct, valid and exclude center")
    delta = xyz[center] - xyz[h]
    scaled_r2 = np.einsum("ij,ij->i", delta, delta) / radius**2
    denominator = 1 + scaled_r2**3
    terms = 1 / denominator
    derivatives = (-6 * scaled_r2**2 / (radius**2 * denominator**2))[:, None] * delta
    gradient = np.zeros_like(xyz)
    gradient[center] = derivatives.sum(axis=0)
    gradient[h] = -derivatives
    return float(terms.sum()), gradient


def coordination_transfer(
    coordinates: FloatArray,
    donor: int,
    acceptor: int,
    hydrogens: Sequence[int],
    radius: float = 1.25,
) -> CoordinationTransfer:
    """Return CN(acceptor,H)−CN(donor,H), independent of hydrogen labels.

    Args:
        coordinates: Finite Cartesian (N,3) array in Å.
        donor: Zero-based donor heavy-atom index (Thr OG1 here).
        acceptor: Zero-based acceptor heavy-atom index (Thr N here).
        hydrogens: All declared real QM hydrogen indices, excluding link atoms.
        radius: Positive switching length in Å.

    Returns:
        Both counts, their difference, and its analytic Cartesian gradient.

    Raises:
        ValueError: If donor equals acceptor or other inputs are invalid.
    """
    if donor == acceptor:
        raise ValueError("Donor and acceptor must differ")
    d_count, d_grad = smooth_coordination(coordinates, donor, hydrogens, radius)
    a_count, a_grad = smooth_coordination(coordinates, acceptor, hydrogens, radius)
    return CoordinationTransfer(d_count, a_count, a_count - d_count, a_grad - d_grad)
