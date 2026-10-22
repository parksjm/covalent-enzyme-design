"""Transfer attached hydrogen coordinates onto updated heavy-atom positions."""

from __future__ import annotations

from typing import Any

import numpy as np


def _minimum_rotation(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Return the deterministic proper rotation mapping one vector to another."""
    source_unit = source / np.linalg.norm(source)
    target_unit = target / np.linalg.norm(target)
    cross = np.cross(source_unit, target_unit)
    sine = float(np.linalg.norm(cross))
    cosine = float(np.clip(np.dot(source_unit, target_unit), -1.0, 1.0))
    if sine < 1e-12:
        if cosine > 0.0:
            return np.eye(3)
        basis = np.eye(3)[int(np.argmin(np.abs(source_unit)))]
        axis = np.cross(source_unit, basis)
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)

    axis = cross / sine
    skew = np.array(
        [[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]]
    )
    return np.eye(3) + sine * skew + (1.0 - cosine) * (skew @ skew)


def _fit_local_rotation(
    old_vectors: np.ndarray, new_vectors: np.ndarray
) -> tuple[np.ndarray, str]:
    """Fit a proper rotation to anchored vectors without centroid alignment."""
    old_lengths = np.linalg.norm(old_vectors, axis=1)
    new_lengths = np.linalg.norm(new_vectors, axis=1)
    if np.any(old_lengths < 1e-12) or np.any(new_lengths < 1e-12):
        raise ValueError("Heavy-atom anchors cannot coincide with their parent")
    old_unit = old_vectors / old_lengths[:, None]
    new_unit = new_vectors / new_lengths[:, None]
    if min(np.linalg.matrix_rank(old_unit), np.linalg.matrix_rank(new_unit)) < 2:
        return _minimum_rotation(old_unit[0], new_unit[0]), "minimum_vector_rotation"

    covariance = new_unit.T @ old_unit
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(left @ right_t)
    rotation = left @ correction @ right_t
    if np.linalg.det(rotation) < 0.0:
        raise ValueError("Local rotation fit produced an improper rotation")
    return rotation, "kabsch"


def transfer_attached_coordinates(
    old_xyz: np.ndarray,
    new_heavy_xyz: dict[int, np.ndarray],
    atomic_numbers: list[int],
    bonds: list[tuple[int, int]],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Move heavy atoms and rotate attached hydrogens with their local frame.

    Args:
        old_xyz: Original coordinates indexed by atom, in Angstroms.
        new_heavy_xyz: Required updated coordinates keyed by heavy-atom index.
        atomic_numbers: Atomic number for each atom in ``old_xyz``.
        bonds: Atom-index pairs describing the molecular connectivity.

    Returns:
        A new coordinate array and a JSON-serializable method/error report.

    Raises:
        ValueError: If coordinates, heavy-atom mappings, or hydrogen parents are
            invalid, or if a hydrogen parent has no heavy-atom anchors.
    """
    coordinates = np.asarray(old_xyz, dtype=float)
    if (
        coordinates.shape != (len(atomic_numbers), 3)
        or not np.isfinite(coordinates).all()
    ):
        raise ValueError("old_xyz must contain one finite three-vector per atom")
    if any(not isinstance(number, int) or number <= 0 for number in atomic_numbers):
        raise ValueError("atomic_numbers must contain positive integers")

    heavy_indices = {i for i, number in enumerate(atomic_numbers) if number != 1}
    if set(new_heavy_xyz) != heavy_indices:
        missing = sorted(heavy_indices - set(new_heavy_xyz))
        extra = sorted(set(new_heavy_xyz) - heavy_indices)
        raise ValueError(
            f"Heavy-atom coordinate mapping mismatch; missing={missing}, extra={extra}"
        )

    new_coordinates = coordinates.copy()
    for index in sorted(heavy_indices):
        position = np.asarray(new_heavy_xyz[index], dtype=float)
        if position.shape != (3,) or not np.isfinite(position).all():
            raise ValueError(
                f"New coordinate for atom {index} must be a finite three-vector"
            )
        new_coordinates[index] = position

    adjacency: list[set[int]] = [set() for _ in atomic_numbers]
    for first, second in bonds:
        if (
            first == second
            or not (0 <= first < len(atomic_numbers))
            or not (0 <= second < len(atomic_numbers))
        ):
            raise ValueError(f"Invalid bond indices: {(first, second)}")
        adjacency[first].add(second)
        adjacency[second].add(first)

    hydrogens_by_parent: dict[int, list[int]] = {}
    for index, number in enumerate(atomic_numbers):
        if number != 1:
            continue
        parents = sorted(
            neighbor for neighbor in adjacency[index] if neighbor in heavy_indices
        )
        if len(parents) != 1:
            raise ValueError(
                f"Hydrogen atom {index} must have exactly one heavy-atom parent"
            )
        hydrogens_by_parent.setdefault(parents[0], []).append(index)

    method_counts = {"kabsch": 0, "minimum_vector_rotation": 0}
    for parent, hydrogens in sorted(hydrogens_by_parent.items()):
        direct = sorted(
            neighbor for neighbor in adjacency[parent] if neighbor in heavy_indices
        )
        anchors = set(direct)
        for neighbor in direct:
            anchors.update(
                candidate
                for candidate in adjacency[neighbor]
                if candidate in heavy_indices and candidate != parent
            )
        anchors.discard(parent)
        ordered_anchors = sorted(anchors)
        if not ordered_anchors:
            raise ValueError(f"Hydrogen parent {parent} has no heavy-atom anchors")
        old_vectors = np.array(
            [coordinates[index] - coordinates[parent] for index in ordered_anchors]
        )
        new_vectors = np.array(
            [
                new_coordinates[index] - new_coordinates[parent]
                for index in ordered_anchors
            ]
        )
        rotation, method = _fit_local_rotation(old_vectors, new_vectors)
        method_counts[method] += len(hydrogens)
        for hydrogen in hydrogens:
            new_coordinates[hydrogen] = new_coordinates[parent] + rotation @ (
                coordinates[hydrogen] - coordinates[parent]
            )

    max_error = 0.0
    for parent, hydrogens in hydrogens_by_parent.items():
        for hydrogen in hydrogens:
            old_length = np.linalg.norm(coordinates[hydrogen] - coordinates[parent])
            new_length = np.linalg.norm(
                new_coordinates[hydrogen] - new_coordinates[parent]
            )
            max_error = max(max_error, abs(float(new_length - old_length)))
    report: dict[str, Any] = {
        "method_counts": method_counts,
        "max_bond_length_error_A": float(max_error),
    }
    return new_coordinates, report
