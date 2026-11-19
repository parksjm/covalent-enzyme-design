"""Parse Gaussian scan windows into energies and signed dihedrals."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from nyl12_covalent.gaussian_checks import final_orientation

SCF_PATTERN = re.compile(r"SCF Done:\s+E\([^)]*\)\s*=\s*([-+\d.Ee]+)")
OPTIMIZATION_END = "Optimization completed."


@dataclass(frozen=True)
class ScanPoint:
    """One completed constrained optimization in a Gaussian scan."""

    index: int
    energy_hartree: float
    angle_degrees: float


def parse_scan_log(
    log_path: Path,
    atom_indices_one_based: list[int],
    expected_points: int,
) -> list[ScanPoint]:
    """Parse completed Gaussian ModRedundant scan points strictly.

    The expected log convention is one ``Optimization completed.`` marker per
    scan window, with the window's final ``SCF Done`` energy and final Cartesian
    orientation printed before that marker. The full log must end in exactly
    one normal Gaussian termination and contain no error termination.

    Args:
        log_path: Gaussian output for one scan direction.
        atom_indices_one_based: Four Gaussian atom indices defining the
            scanned dihedral, in input order.
        expected_points: Number of constrained geometries expected in the log.

    Returns:
        Ordered scan points with final energy and optimized dihedral.

    Raises:
        FileNotFoundError: If the Gaussian output is absent.
        ValueError: If output is incomplete, malformed, or has a different
            number of windows than declared.
    """
    if len(atom_indices_one_based) != 4 or min(atom_indices_one_based) < 1:
        raise ValueError("A scan dihedral requires four positive one-based indices")
    if expected_points < 1:
        raise ValueError("expected_points must be positive")
    text = log_path.read_text(encoding="utf-8", errors="replace")
    if text.count("Normal termination of Gaussian 16") != 1:
        raise ValueError(f"Expected one normal Gaussian termination: {log_path}")
    if "Error termination" in text:
        raise ValueError(f"Gaussian reported an error termination: {log_path}")
    ends = list(re.finditer(re.escape(OPTIMIZATION_END), text))
    if len(ends) != expected_points:
        raise ValueError(
            f"Expected {expected_points} completed scan windows, found {len(ends)}"
        )
    points: list[ScanPoint] = []
    segment_start = 0
    for index, end in enumerate(ends, start=1):
        segment = text[segment_start : end.start()]
        energies = SCF_PATTERN.findall(segment)
        if not energies:
            raise ValueError(f"No SCF energy for scan window {index}: {log_path}")
        elements, coordinates = final_orientation(segment)
        if max(atom_indices_one_based) > len(elements):
            raise ValueError(f"Dihedral index exceeds atom count in {log_path}")
        atom_indices = [value - 1 for value in atom_indices_one_based]
        angle = dihedral_degrees(coordinates[atom_indices])
        points.append(ScanPoint(index, float(energies[-1]), angle))
        segment_start = end.end()
    return points


def dihedral_degrees(coordinates: np.ndarray) -> float:
    """Calculate a signed dihedral in degrees for four Cartesian points.

    Args:
        coordinates: Array of shape ``(4, 3)`` in Angstrom.

    Returns:
        Signed angle in the interval ``[-180, 180]``.

    Raises:
        ValueError: If points are malformed or define a degenerate plane.
    """
    points = np.asarray(coordinates, dtype=float)
    if points.shape != (4, 3) or not np.isfinite(points).all():
        raise ValueError("Dihedral coordinates must be a finite 4 by 3 array")
    first, second, third, fourth = points
    axis = third - second
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-10:
        raise ValueError("Dihedral central bond has zero length")
    axis /= axis_norm
    vector_a = first - second
    vector_b = fourth - third
    vector_a -= np.dot(vector_a, axis) * axis
    vector_b -= np.dot(vector_b, axis) * axis
    if np.linalg.norm(vector_a) < 1e-10 or np.linalg.norm(vector_b) < 1e-10:
        raise ValueError("Dihedral atoms are collinear")
    x_value = float(np.dot(vector_a, vector_b))
    y_value = float(np.dot(np.cross(axis, vector_a), vector_b))
    return math.degrees(math.atan2(y_value, x_value))
