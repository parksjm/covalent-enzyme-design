"""Tests of gate attribution and deterministic representative selection."""

import numpy as np
import pytest
from nyl12_covalent.native_thr_md_gate_audit import representatives, summarize, values


def row(
    o: float,
    n: float,
    p: float,
    c: float,
    site: str = "AE",
    time: int = 1,
    tier: str = "A",
) -> dict:
    """Make a saved-frame-like scalar row."""
    return {
        "og1_shift": o,
        "nterm_shift": n,
        "psi_delta": p,
        "chi1_delta": c,
        "scaffold_rmsd": 0.2,
        "site": site,
        "time_ps": time,
        "tier": tier,
    }


def test_attribution_and_diagnostics() -> None:
    data = [row(0.5, 1, 20, 10), row(1, 1, 20, 10), row(1, 2, 50, 10)]
    s = summarize(data)
    assert s["n"] == 3 and s["fails"]["og1_shift"] == 2
    assert s["exclusive_fails"]["og1_shift"] == 1
    assert s["leave_one_out_passes"]["og1_shift"] == 2
    assert s["failure_combinations"]["og1_shift+nterm_shift+psi_delta"] == 1
    assert s["one_at_a_time_diagnostic"]["og1_shift"][1]["passes"] == 2
    assert s["joint_diagnostic"][1]["passes"] == 3


def test_representatives_and_invalid_inputs() -> None:
    data = [
        row(0.4, 0.8, 10, 5, time=1),
        row(0.5, 0.7, 9, 4, time=2),
        row(0.2, 0.5, 5, 2, site="BF", time=3),
    ]
    s = representatives(data)
    assert s["closest_native"]["site"] == "BF"
    assert s["AE_median_medoid"]["site"] == "AE"
    with pytest.raises(ValueError):
        summarize([])
    with pytest.raises(ValueError):
        values([row(np.nan, 1, 1, 1)])
    with pytest.raises(ValueError):
        representatives([row(1, 1, 1, 1, site="BF")])
