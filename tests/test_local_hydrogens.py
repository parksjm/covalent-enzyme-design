"""Tests for local attached-hydrogen coordinate transfer."""

import numpy as np
import pytest

from nyl12_covalent.local_hydrogens import (
    _fit_local_rotation,
    transfer_attached_coordinates,
)


def _coordinates() -> tuple[
    np.ndarray, dict[int, np.ndarray], list[int], list[tuple[int, int]]
]:
    old = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ]
    )
    new = {
        0: np.array([0.3, -0.2, 0.1]),
        1: np.array([1.3, -0.2, 0.1]),
        2: np.array([0.3, 0.8, 0.1]),
        3: np.array([0.3, -0.2, 1.1]),
    }
    return old, new, [6, 6, 6, 6, 1], [(0, 1), (0, 2), (0, 3), (0, 4)]


def test_rank_two_frame_is_rigid_transform_equivariant() -> None:
    old, new, numbers, bonds = _coordinates()
    moved, report = transfer_attached_coordinates(old, new, numbers, bonds)
    angle = 0.7
    rigid = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    shift = np.array([-1.0, 2.0, 0.4])
    old_rigid = old @ rigid.T + shift
    new_rigid = {index: position @ rigid.T + shift for index, position in new.items()}
    moved_rigid, _ = transfer_attached_coordinates(old_rigid, new_rigid, numbers, bonds)
    assert moved_rigid == pytest.approx(moved @ rigid.T + shift)
    assert report["method_counts"]["kabsch"] == 1
    assert report["max_bond_length_error_A"] < 1e-12


def test_local_heavy_atom_movement_rotates_hydrogen() -> None:
    old, new, numbers, bonds = _coordinates()
    new[2] = np.array([0.3, 0.0, 1.1])
    moved, _ = transfer_attached_coordinates(old, new, numbers, bonds)
    assert moved[4] != pytest.approx(new[0] + (old[4] - old[0]))
    assert np.linalg.norm(moved[4] - new[0]) == pytest.approx(
        np.linalg.norm(old[4] - old[0])
    )


def test_reflected_anchor_frame_uses_proper_rotation() -> None:
    old, new, numbers, bonds = _coordinates()
    new[3] = np.array([0.0, 0.0, -1.0])
    moved, _ = transfer_attached_coordinates(old, new, numbers, bonds)
    assert np.linalg.norm(moved[4] - new[0]) == pytest.approx(1.0)
    rotation, _ = _fit_local_rotation(
        old[1:4] - old[0], np.array([new[i] - new[0] for i in (1, 2, 3)])
    )
    assert np.linalg.det(rotation) == pytest.approx(1.0)
    assert rotation.T @ rotation == pytest.approx(np.eye(3))


@pytest.mark.parametrize(
    "target", [np.array([-1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])]
)
def test_single_anchor_minimum_rotation_including_antiparallel(
    target: np.ndarray,
) -> None:
    old = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    new = {0: np.zeros(3), 1: target}
    moved, report = transfer_attached_coordinates(old, new, [6, 6, 1], [(0, 1), (0, 2)])
    assert moved[2] == pytest.approx(
        _rotation_for_test(target) @ np.array([0.0, 1.0, 0.0])
    )
    assert report["method_counts"]["minimum_vector_rotation"] == 1


def _rotation_for_test(target: np.ndarray) -> np.ndarray:
    if target[0] < 0:
        return np.diag([-1.0, -1.0, 1.0])
    return np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])


@pytest.mark.parametrize(
    ("numbers", "bonds", "mapping", "message"),
    [
        ([6, 1], [(0, 1)], {}, "mapping mismatch"),
        ([6, 1], [(1, 1)], {0: np.zeros(3)}, "Invalid bond"),
        ([6, 1], [], {0: np.zeros(3)}, "exactly one heavy-atom parent"),
        ([6, 1], [(0, 1)], {0: np.zeros(3)}, "no heavy-atom anchors"),
    ],
)
def test_invalid_inputs_fail_explicitly(
    numbers: list[int],
    bonds: list[tuple[int, int]],
    mapping: dict[int, np.ndarray],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        transfer_attached_coordinates(
            np.zeros((len(numbers), 3)), mapping, numbers, bonds
        )
