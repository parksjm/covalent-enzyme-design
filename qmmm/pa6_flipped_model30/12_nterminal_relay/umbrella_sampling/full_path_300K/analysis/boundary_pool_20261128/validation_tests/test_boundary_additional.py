"""Actual baseline diagnostic validation and unavailable-cohort fail-closed gates."""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pytest
import boundary_collect as module
from diagnostics import summarize
import watch  # noqa: F401 - expose orchestration separately in measured coverage


def test_actual_frozen_baseline_geometry_order_and_support() -> None:
    """Replay one actual1300-frame history with all direct chemistry identities."""
    base = Path(module.registration()["baseline"])
    sources = json.loads((base / "sources.json").read_text())
    weights = np.load(base / "weights.npz")
    state = sources["states"][0]
    count = int(weights["counts"][0])
    data = dict(
        sources,
        states=[state],
        proofs=[p for p in sources["proofs"] if p["state_id"] == state["state_id"]],
        q=weights["q"][:count],
    )
    result = summarize(data, weights["logweights"][:count])
    assert result["direct_CV_identities_verified"] == 3
    assert result["exact_weight_order"]
    assert result["histories"][0]["state_id"] == state["state_id"]
    assert len(result["histories"][0]["distances"]) == 5
    assert result["support"]
    assert all(p["weight_Kish_ESS"] > 0 for p in result["support"])


def test_collect_rejects_missing_actual_boundary_stages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A vanished/incomplete accepted prefix fails before combining any histories."""
    cfg = module.registration()
    source = json.loads((Path(cfg["baseline"]) / "sources.json").read_text())
    monkeypatch.setattr(
        module, "collect", lambda selected, output: dict(config=source["config"])
    )
    with pytest.raises(ValueError, match="prefix disappeared"):
        module.collect_boundary(30, tmp_path)
