"""Differentiate frozen dihedrals and construct Cartesian rigid motions."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

FloatArray = NDArray[np.float64]
Torsion = tuple[int, int, int, int]


def _angle(points: FloatArray) -> float:
    """Measure a nondegenerate four-atom torsion in radians."""
    edges = np.diff(points, axis=0)
    lengths = np.linalg.norm(edges, axis=1)
    if np.any(lengths < 1e-10):
        raise ValueError("Coincident atoms in a torsion constraint")
    unit = edges / lengths[:, None]
    if (
        min(
            np.linalg.norm(np.cross(unit[0], unit[1])),
            np.linalg.norm(np.cross(unit[1], unit[2])),
        )
        < 1e-8
    ):
        raise ValueError("Collinear atoms in a torsion constraint")
    conformer = Chem.Conformer(4)
    for index, point in enumerate(points):
        conformer.SetAtomPosition(index, tuple(float(value) for value in point))
    angle = float(rdMolTransforms.GetDihedralRad(conformer, 0, 1, 2, 3))
    if not np.isfinite(angle):
        raise ValueError("Nonfinite torsion angle")
    return angle


def dihedral_constraint_derivatives(
    coordinates_bohr: FloatArray,
    torsions: Sequence[Torsion],
    *,
    step_bohr: float,
) -> tuple[FloatArray, FloatArray]:
    """Central-difference the local wrapped torsion constraints.

    Only each torsion's twelve Cartesian coordinates are displaced. Angles
    are wrapped relative to the undisplaced value to avoid a +/-pi branch cut.
    The caller must compare multiple step sizes before accepting derivatives.

    Args:
        coordinates_bohr: Cartesian coordinates, shape (N, 3), in bohr.
        torsions: Unique ordered four-tuples of zero-based atom indexes.
        step_bohr: Explicit finite-difference Cartesian displacement in bohr.

    Returns:
        Jacobian (m, 3N) in radian/bohr and constraint Hessians (m, 3N, 3N)
        in radian/bohr^2. Neither contains energy derivatives.

    Raises:
        ValueError: For invalid coordinates, torsions, step or degeneracies.
    """
    xyz = np.asarray(coordinates_bohr, dtype=float)
    if (
        xyz.ndim != 2
        or xyz.shape[1] != 3
        or xyz.shape[0] < 4
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Expected finite (N,3) coordinates with N >= 4")
    if not np.isfinite(step_bohr) or step_bohr <= 0:
        raise ValueError("Finite-difference step must be finite and positive")
    if not torsions or len(set(torsions)) != len(torsions):
        raise ValueError("Expected nonempty, unique torsion definitions")
    jacobian = np.zeros((len(torsions), xyz.size))
    hessians = np.zeros((len(torsions), xyz.size, xyz.size))
    for row, torsion in enumerate(torsions):
        if len(torsion) != 4 or any(
            not isinstance(i, (int, np.integer)) for i in torsion
        ):
            raise ValueError("Torsion requires four integer atom indexes")
        if len(set(torsion)) != 4 or min(torsion) < 0 or max(torsion) >= len(xyz):
            raise ValueError("Torsion atom indexes must be distinct and in range")
        points = xyz[list(torsion)].copy()
        origin = _angle(points)

        def relative(displacement: FloatArray) -> float:
            delta = _angle(points + displacement.reshape((4, 3))) - origin
            return float(np.arctan2(np.sin(delta), np.cos(delta)))

        global_indexes = np.asarray(
            [3 * atom + axis for atom in torsion for axis in range(3)]
        )
        for first in range(12):
            offset = np.zeros(12)
            offset[first] = step_bohr
            plus, minus = relative(offset), relative(-offset)
            first_global = global_indexes[first]
            jacobian[row, first_global] = (plus - minus) / (2 * step_bohr)
            hessians[row, first_global, first_global] = (plus + minus) / step_bohr**2
            for second in range(first):
                other = np.zeros(12)
                other[second] = step_bohr
                value = (
                    relative(offset + other)
                    - relative(offset - other)
                    - relative(-offset + other)
                    + relative(-offset - other)
                ) / (4 * step_bohr**2)
                second_global = global_indexes[second]
                hessians[row, first_global, second_global] = value
                hessians[row, second_global, first_global] = value
    return jacobian, hessians


def rigid_cartesian_directions(coordinates_bohr: FloatArray) -> FloatArray:
    """Construct independent translations and rotations of a nonlinear molecule.

    Args:
        coordinates_bohr: Finite Cartesian coordinates with shape (N, 3).

    Returns:
        Six unnormalized Cartesian rigid-motion columns, shape (3N, 6).
        The arithmetic centroid fixes an origin; masses are not used.

    Raises:
        ValueError: If coordinates are invalid or the molecule is linear.
    """
    xyz = np.asarray(coordinates_bohr, dtype=float)
    if (
        xyz.ndim != 2
        or xyz.shape[1] != 3
        or xyz.shape[0] < 3
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Expected finite nonlinear-molecule (N,3) coordinates")
    centered = xyz - xyz.mean(axis=0)
    translations = [np.tile(axis, (len(xyz), 1)).ravel() for axis in np.eye(3)]
    rotations = [np.cross(axis, centered).ravel() for axis in np.eye(3)]
    result = np.column_stack(translations + rotations)
    singular = np.linalg.svd(result, compute_uv=False)
    if np.count_nonzero(singular > singular[0] * 1e-10) != 6:
        raise ValueError("Six rigid motions require a nonlinear molecule")
    return result
