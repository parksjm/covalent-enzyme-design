"""Test actual MBAR refits, reproducibility, and lost-support accounting."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.relay_pmf_block_bootstrap import block_bootstrap


def arguments() -> dict:
    """Make two identical umbrellas, four blocks each, with one rare bin."""
    s = np.array([0.5] + [1.5] * 39 + [2.5] * 40, float)
    return dict(
        reduced=np.zeros((2, 80)),
        counts=np.array([40, 40]),
        s=s,
        edges=np.array([0, 1, 2, 3, 4]),
        temperature=300.0,
        reference_s=1.5,
        block_lengths=[[10] * 4, [10] * 4],
        replicates=50,
        seed=16,
        workers=1,
    )


def test_real_refits_are_reproducible_and_do_not_hide_lost_support() -> None:
    first, second = block_bootstrap(**arguments()), block_bootstrap(**arguments())
    assert first == second
    assert 0 < first["support_fraction"][0] < 1
    assert first["support_fraction"][1:3] == [1.0, 1.0]
    assert first["support_fraction"][3] == 0
    assert first["lower_95_kcal_mol"][3] is None
    assert first["lower_95_kcal_mol"][1] == pytest.approx(0)
    assert sum(row[0] is None for row in first["replicates_kcal_mol"]) > 0


def test_process_pool_matches_serial_seeded_results() -> None:
    serial = block_bootstrap(**arguments())
    parallel = block_bootstrap(**(arguments() | dict(workers=2)))
    for key in ("replicates_kcal_mol", "support_fraction", "lower_95_kcal_mol"):
        assert serial[key] == parallel[key]


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (dict(replicates=3), "replicates"),
        (dict(workers=9), "workers"),
        (dict(counts=np.array([10, 10])), "matching"),
        (dict(temperature=-1), "temperature"),
        (dict(edges=np.array([0, 0, 4])), "edges"),
        (dict(reference_s=10), "Reference"),
        (dict(block_lengths=[[40], [40]]), "two nonempty"),
    ],
)
def test_invalid_bootstrap_inputs_fail(change: dict, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        block_bootstrap(**(arguments() | change))
