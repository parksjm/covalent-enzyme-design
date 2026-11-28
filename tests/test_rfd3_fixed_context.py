"""Enforce the narrow immutable terminal-oxygen normalization contract."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.rfd3_fixed_context import reconcile_terminal_oxygens


def fixture_atoms() -> tuple[list, np.ndarray, list, np.ndarray]:
    """Create eight distinct rigid terminal contexts with one omitted OXT each."""
    keys = [
        (c, 232 if c in "ABCD" else 96, n)
        for c in "ABCDEFGH"
        for n in ("C", "O", "OXT")
    ]
    rng = np.random.default_rng(7)
    source = rng.normal(size=(24, 3)) * 10
    native = [k for k in keys if k[2] != "OXT"]
    xyz = source[[i for i, k in enumerate(keys) if k[2] != "OXT"]] + 5
    return keys, source, native, xyz


def test_exact_source_restoration() -> None:
    """Only omitted source OXT positions are copied after proper rigid alignment."""
    keys, source, native, xyz = fixture_atoms()
    restored, receipt = reconcile_terminal_oxygens(keys, source, native, xyz, set())
    assert np.allclose(restored, source)
    assert len(receipt["restored_source_atoms"]) == 8
    assert receipt["raw_native_coordinates_modified"] is False
    assert np.allclose(receipt["raw_to_source_transform"]["translation"], [-5, -5, -5])


def test_reject_unexpected_missing_mobile_duplicate_and_moved_anchors() -> None:
    """Never reconstruct unexpected atoms or a moving terminal context."""
    keys, source, native, xyz = fixture_atoms()
    with pytest.raises(ValueError, match="exactly eight"):
        reconcile_terminal_oxygens(keys, source, native[:-1], xyz[:-1], set())
    with pytest.raises(ValueError, match="Ambiguous"):
        reconcile_terminal_oxygens(
            keys, source, [*native, native[0]], np.vstack([xyz, xyz[0]]), set()
        )
    with pytest.raises(ValueError, match="mobile"):
        reconcile_terminal_oxygens(keys, source, native, xyz, {"A232"})
    changed = xyz.copy()
    changed[0] += 1
    with pytest.raises(ValueError, match="anchor moved"):
        reconcile_terminal_oxygens(keys, source, native, changed, set())
