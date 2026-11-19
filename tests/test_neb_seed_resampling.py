"""Focused tests for trajectory-based native NEB seed resampling."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb_prepare.py"
)
SPEC = importlib.util.spec_from_file_location("new_neb_prepare", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_resample_path_keeps_endpoints_and_uniform_arc_positions() -> None:
    """Resampling retains the endpoints and spaces a straight path evenly."""
    frames = [
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[1.0, 0.0, 0.0]]),
        np.array([[3.0, 0.0, 0.0]]),
    ]

    metric = [np.array([frame[0, 0]]) for frame in frames]
    result = MODULE.resample_path(frames, metric, intervals=3)

    assert len(result) == 4
    np.testing.assert_array_equal(result[0], frames[0])
    np.testing.assert_array_equal(result[-1], frames[-1])
    np.testing.assert_allclose([frame[0, 0] for frame in result], [0.0, 1.0, 2.0, 3.0])


def test_resample_path_rejects_degenerate_or_malformed_input() -> None:
    """Invalid paths fail explicitly instead of yielding duplicate images."""
    with pytest.raises(ValueError, match="zero arc length"):
        MODULE.resample_path(
            [np.zeros((1, 3)), np.zeros((1, 3))],
            [np.zeros(1), np.zeros(1)],
            intervals=2,
        )
    with pytest.raises(ValueError, match="matching finite"):
        MODULE.resample_path(
            [np.zeros((1, 3)), np.ones((2, 3))],
            [np.zeros(1), np.ones(1)],
            intervals=2,
        )


def test_build_writes_regular_neb_inputs_without_launching(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The full real-source preflight writes reproducible inputs only."""
    monkeypatch.setattr(MODULE, "OUTPUT", tmp_path / "new_neb")

    config = MODULE.build()

    assert config["seed_structure_count"] == 13
    assert config["seed_ts_index_zero_based"] == 6
    assert config["ci_enabled"] is False
    assert config["launched"] is False
    assert max(config["maximum_adjacent_reactive_distance_delta_A"].values()) < 0.36
    assert len(list((tmp_path / "new_neb/inputs/images").glob("image_*.ncrst"))) == 13
    assert (
        "neb_climb_test=0.0"
        in (tmp_path / "new_neb/calculation/minimize.mdin").read_text()
    )
