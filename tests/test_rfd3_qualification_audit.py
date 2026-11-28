"""Check source-mapped RFD3 acceptance independently of native inference."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.rfd3_qualification_audit import align_on_fixed, coordinate_gates


def test_proper_alignment_and_invalid_inputs() -> None:
    """Rigid placement differences vanish without permitting a mirror transform."""
    source = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 1, 0], [0.0, 0, 1]])
    fixed = np.ones(4, dtype=bool)
    assert np.allclose(align_on_fixed(source, source + 10, fixed), source)
    with pytest.raises(ValueError, match="shape"):
        align_on_fixed(source, source[:2], fixed)
    with pytest.raises(ValueError, match="Invalid"):
        align_on_fixed(source, source * np.nan, fixed)


def test_gates_distinguish_no_movement_clash_closure_and_chirality() -> None:
    """Independent gates reject stationary, distorted, mirrored and clashing output."""
    atoms = [{"chain": "A", "residue": 1, "name": n} for n in ("N", "CA", "C", "CB")]
    atoms += [{"chain": "A", "residue": 2, "name": "N"}]
    atoms += [{"chain": "B", "residue": 1, "name": n} for n in ("N", "CA", "C", "CB")]
    source = np.array(
        [
            [0.0, 0, 0],
            [1.0, 0, 0],
            [1.0, 1, 0],
            [1.0, 0, 1],
            [2.3, 1, 0],
            [20.0, 20, 20],
            [21.0, 20, 20],
            [21.0, 21, 20],
            [21.0, 20, 21],
        ]
    )
    edges = [[0, 1], [1, 2], [1, 3], [2, 4], [5, 6], [6, 7], [6, 8]]
    thresholds = {
        "fixed_max_A": 0.001,
        "mobile_CA_rmsd_min_A": 0.01,
        "mobile_CA_rmsd_max_A": 3.0,
        "peptide_min_A": 1.15,
        "peptide_max_A": 1.65,
        "clash_min_A": 0.3,
    }

    def audit(output: np.ndarray) -> dict:
        return coordinate_gates(
            source, output, atoms, edges, ["A1"], [0, 8], thresholds
        )

    moved = source.copy()
    moved[1, 0] += 0.1
    assert audit(moved)["accepted"]
    assert not audit(source)["accepted"]
    mirrored = moved.copy()
    mirrored[3, 2] = -1
    assert audit(mirrored)["stereochemical_inversions"] == ["A1:CA"]
    broken = moved.copy()
    broken[2] = [1.0, 10, 0]
    assert not audit(broken)["peptide_closure_pass"]
    clashing = moved.copy()
    clashing[0] = source[5]
    assert audit(clashing)["severe_heavy_clashes"]
    badfixed = moved.copy()
    badfixed[4, 0] += 1
    assert not audit(badfixed)["accepted"]
    with pytest.raises(ValueError, match="Empty"):
        coordinate_gates(source, source, atoms, edges, [], [0, 8], thresholds)
    with pytest.raises(ValueError, match="counts"):
        coordinate_gates(source, source, atoms[:2], edges, ["A1"], [0, 8], thresholds)
    duplicate = [*atoms[:-1], atoms[0]]
    with pytest.raises(ValueError, match="Ambiguous"):
        coordinate_gates(source, source, duplicate, edges, ["A1"], [0, 8], thresholds)
