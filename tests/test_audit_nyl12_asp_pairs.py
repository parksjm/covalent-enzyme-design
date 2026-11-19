"""Scientific checks for symmetric Asp oxygen distances."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "asp_audit", Path(__file__).parents[1] / "scripts/audit_nyl12_asp_pairs.py"
)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_pairs_and_label_invariance() -> None:
    xyz = np.array([[[0, 0, 0], [0, 1, 0], [3, 0, 0], [3, 2, 0]]])
    measured = AUDIT.pair_distances(xyz, None)
    np.testing.assert_allclose(measured, [[3, np.sqrt(13), np.sqrt(10), np.sqrt(10)]])
    swapped = AUDIT.pair_distances(xyz[:, [1, 0, 3, 2]], None)
    np.testing.assert_allclose(measured.min(axis=1), swapped.min(axis=1))


def test_periodic_crossing() -> None:
    xyz = np.array([[[0.5, 0, 0], [0.5, 1, 0], [9.5, 0, 0], [9.5, 1, 0]]])
    measured = AUDIT.pair_distances(xyz, np.array([10, 10, 10, 90, 90, 90]))
    np.testing.assert_allclose(measured, [[1, np.sqrt(2), np.sqrt(2), 1]])


def test_truncated_octahedron_translation() -> None:
    angle = np.deg2rad(109.471219)
    cosine, sine = np.cos(angle), np.sin(angle)
    cy = (cosine - cosine**2) / sine
    cell = 10 * np.array([[1, 0, 0], [cosine, sine, 0],
                          [cosine, cy, np.sqrt(1 - cosine**2 - cy**2)]])
    xyz = np.array([[[0, 0, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0]]], dtype=float)
    expected = AUDIT.pair_distances(xyz, None)
    xyz[:, 2:] += 2 * cell[0] - 3 * cell[1] + cell[2]
    measured = AUDIT.pair_distances(xyz, np.array([10, 10, 10, *([109.471219] * 3)]))
    np.testing.assert_allclose(measured, expected, atol=1e-6)


def test_invalid_box_fails() -> None:
    with pytest.raises(ValueError):
        AUDIT.pair_distances(np.zeros((1, 4, 3)), np.array([10, 10, 0, 90, 90, 90]))


def test_summary_minimum_per_frame() -> None:
    values = np.array([[[2, 4, 5, 6]], [[6, 5, 4, 2]]])
    summary = AUDIT.write_summary("replica", values)
    assert summary[-1]["pair"] == "minimum"
    assert summary[-1]["mean_A"] == 2
    assert summary[-1]["fraction_lt_3A"] == 1
