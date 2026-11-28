"""Check extended physical-source selection and one-time production discard."""

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from nyl12_covalent.relay_fixed_burnin import fixed_window_data
from nyl12_covalent.relay_mbar_progression import select_horizon

ROOT = Path(__file__).resolve().parents[1]
BASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K"
)
SIDE = BASE / "analysis/main100_continuation_20261127"


@pytest.mark.parametrize("horizon", [65, 70, 75, 80, 85, 90, 95, 100])
def test_extended_native_blocks_and_single_discard(
    tmp_path: Path, horizon: int
) -> None:
    """Retain block13onwards and discard the first5ps only once per lineage."""
    paths, blocks = [], []
    for i in range(horizon // 5):
        p = tmp_path / "one" / f"{i + 3:02d}_production_5ps/analysis_acceptance.json"
        p.parent.mkdir(parents=True)
        p.write_text(
            json.dumps(dict(window="one", accepted=True, sampled_duration_ps=5))
        )
        paths.append(p)
        times = 6 + i * 5 + np.arange(1, 101) * 0.05
        blocks.append(
            dict(
                window="one",
                sampled_duration_ps=5,
                values=np.zeros((100, 3)),
                original_rows=[dict(time_ps=t) for t in times],
            )
        )
    selected, durations = select_horizon(paths, [dict(label="one")], horizon)
    assert len(selected) == horizon // 5
    assert durations == dict(one=horizon)
    _, counts, selection = fixed_window_data([dict(label="one")], blocks, 5, horizon)
    assert counts.tolist() == [(horizon - 5) * 20]
    assert selection[0]["selected_start_ps"] == pytest.approx(5.05)
    assert select_horizon(paths[:-1], [dict(label="one")], horizon) is None


@pytest.mark.parametrize("horizon", [65, 100])
def test_copied_checkpoint_guard_accepts_future_incomplete(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, horizon: int
) -> None:
    """Future accepted blocks are recognized without claiming incomplete sources."""
    sys.path.insert(0, str(SIDE / "source_snapshot/scripts"))
    spec = importlib.util.spec_from_file_location(
        "future_checkpoint_sources", SIDE / "checkpoint_sources.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    old = tmp_path / "analysis/fixed_burnin_20261126/checkpoint_35ps_refresh_20261126"
    old.mkdir(parents=True)
    (old / "analysis_request.json").write_text(
        json.dumps(dict(windows=[dict(label="one")]))
    )
    monkeypatch.setattr(module, "discover", lambda *_: ([], False))
    assert module.bind_checkpoint(tmp_path, tmp_path / "output", horizon) is None
    with pytest.raises(ValueError, match="authorized"):
        module.bind_checkpoint(tmp_path, tmp_path / "output", 60)
