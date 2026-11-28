"""Focused inventory, native-prefix and fixed-discard regressions."""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pytest
import boundary_collect as module


def synthetic(n: int, start: int, frames: int) -> dict:
    """Build independent histories with controlled duplicate bias identities."""
    return dict(
        states=[
            dict(state_id=str(i), bias=dict(center=i % 95))
            for i in range(start, start + n)
        ],
        proofs=[
            dict(source_sha256={f"{i}.csv": str(i)}) for i in range(start, start + n)
        ],
        counts=np.repeat(frames, n),
        q=np.zeros((frames * n, 3)),
    )


def test_combine_exact_inventory_and_once_discard() -> None:
    """Keep duplicate-bias starts distinct and retain exactly5ps after10ps."""
    result = module.combine(synthetic(100, 0, 10), synthetic(6, 100, 100), 10)
    assert len(result["states"]) == 106
    assert result["counts"][-6:].tolist() == [100] * 6
    assert len({json.dumps(s["bias"]) for s in result["states"]}) == 95


def test_wrong_discard_fails() -> None:
    """Reject an incorrect twice-discarded or undiscarded history count."""
    with pytest.raises(ValueError, match="discard"):
        module.combine(synthetic(100, 0, 10), synthetic(6, 100, 200), 10)


def test_duplicated_block_fails() -> None:
    """Reject recycled physical CSV identities across cohorts."""
    new = synthetic(6, 100, 100)
    new["proofs"][0]["source_sha256"] = {"new.csv": "0"}
    with pytest.raises(ValueError, match="Duplicate"):
        module.combine(synthetic(100, 0, 10), new, 10)


def test_ready_requires_all_consecutive_production_stages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Equilibration and five complete histories cannot stand in for a sixth."""
    prepared = tmp_path / "prepared.json"
    prepared.write_text(
        json.dumps(dict(windows=[dict(label=str(i)) for i in range(6)]))
    )
    config = dict(
        boundary_prepared=str(prepared),
        boundary_prepared_sha256=module.digest(prepared),
        minimum_new_production_checkpoints_ps=[10, 15],
    )
    monkeypatch.setattr(module, "registration", lambda: config)
    for i in range(6):
        for number in (1, 2, 3, 4):
            stage = tmp_path / str(i) / f"{number:02d}_production_5ps"
            stage.mkdir(parents=True)
            (stage / "analysis_acceptance.json").write_text("{}")
    missing = tmp_path / "5/04_production_5ps/analysis_acceptance.json"
    missing.unlink()
    assert not module.ready(10)
    missing.write_text("{}")
    assert module.ready(10)
    assert not module.ready(15)
    with pytest.raises(ValueError, match="Unauthorized"):
        module.ready(5)
    prepared.write_text("{}")
    with pytest.raises(ValueError, match="identity"):
        module.ready(10)
