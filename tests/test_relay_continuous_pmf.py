"""Verify temporal concatenation, single burn-in, and source-bound MBAR refits."""

from __future__ import annotations

from copy import deepcopy

import numpy as np
import pytest

from nyl12_covalent import relay_continuous_pmf as module
from nyl12_covalent.relay_pmf_block_bootstrap import block_bootstrap


def fixture() -> tuple[dict, list[dict], dict]:
    """Create two windows with continuous segment timestamps and exact biases."""
    rng = np.random.default_rng(45)
    windows, blocks = [], []
    for index, center in enumerate([-0.3, 0.3]):
        label = f"window{index}"
        windows.append(
            dict(
                label=label,
                bias=dict(
                    feature_names=["q_acyl", "q_act"],
                    centers_A=[center, 0],
                    force_constants=[1, 1],
                    outer_extent_A=10,
                ),
            )
        )
        for segment in range(4):
            values = rng.normal([center / 2, 0, 0], [0.3, 0.3, 0.4], (100, 3))
            times = 0.05 * (np.arange(100) + segment * 100 + 1)
            blocks.append(
                dict(
                    window=label,
                    values=values,
                    original_rows=[dict(time_ps=str(value)) for value in times],
                    sampled_duration_ps=5,
                    receipt_sha256=f"{len(blocks):064x}",
                    source_sha256={"source.csv": "a" * 64},
                )
            )
    config = dict(
        windows=windows,
        temperature_K=300,
        path_3cv_A=[[-3, 0, 0], [3, 0, 0]],
        s_bin_edges=np.linspace(-0.1, 1.1, 25).tolist(),
        es_search_bounds=[0.35, 0.55],
        acyl_search_bounds=[0.55, 0.7],
        neighbor_pairs=[[0, 1]],
        physical_hamiltonian_sha256="b" * 64,
    )
    parent = dict(
        accepted_sources=[
            dict(
                receipt_sha256=block["receipt_sha256"],
                source_sha256=block["source_sha256"],
            )
            for block in blocks
        ],
        temperature_K=300,
        checkpoint=dict(plot_common_es_reference_s=0.475),
        path_3cv_sha256="c" * 64,
        gates=dict(neighbor_histogram_overlap=dict(minimum=0.2)),
        analysis_revision=dict(old=True),
    )
    return config, blocks, parent


def test_segments_sorted_then_selected_once_per_window(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config, blocks, _ = fixture()
    calls = []

    def select(values: np.ndarray) -> tuple[np.ndarray, dict]:
        calls.append(values.copy())
        return np.arange(150, 400, 5), dict(
            raw_frames=400,
            discarded_frames=150,
            statistical_inefficiency=5,
            selected_frames=50,
        )

    monkeypatch.setattr(module, "select_uncorrelated", select)
    samples, counts, summaries = module.continuous_selections(
        config["windows"], blocks[::-1]
    )
    assert len(calls) == 2
    np.testing.assert_equal(
        calls[0], np.concatenate([part["values"] for part in blocks[:4]])
    )
    assert counts.tolist() == [50, 50]
    assert samples.shape == (100, 3)
    assert summaries[0]["burn_in_ps"] == pytest.approx(7.5)


@pytest.mark.parametrize("kind", ["gap", "duplicate", "unknown", "missing"])
def test_temporal_gaps_and_unknown_or_missing_windows_fail(kind: str) -> None:
    config, blocks, _ = fixture()
    if kind == "gap":
        blocks[1]["original_rows"][0]["time_ps"] = "5.2"
    elif kind == "duplicate":
        blocks.append(deepcopy(blocks[0]))
    elif kind == "unknown":
        blocks[0]["window"] = "unknown"
    else:
        blocks = blocks[:4]
    with pytest.raises(ValueError):
        module.continuous_selections(config["windows"], blocks)


def test_actual_continuous_fit_preserves_source_and_has_explicit_intervals() -> None:
    config, blocks, parent = fixture()
    original = deepcopy(parent)
    report, uncertainty = module.analyze_continuous(
        config, blocks, parent, replicates=50, workers=1
    )
    assert parent == original
    assert report["selection_policy"] == "continuous_window"
    assert len(report["window_selection"]) == 2
    assert len(report["blocks"]) == 8
    assert report["converged"] is False
    assert "analysis_revision" not in report
    assert (
        uncertainty["window_bootstrap"]["resampling_unit"]
        == "decorrelated_configuration"
    )
    assert uncertainty["reference_s"] == 0.475
    assert np.array(report["overlap"]).shape == (2, 2)
    assert sum(report["N_k"]) == sum(uncertainty["selected_frames_per_bin"])
    bad = deepcopy(parent)
    bad["accepted_sources"][0]["receipt_sha256"] = "bad"
    with pytest.raises(ValueError, match="accepted parent"):
        module.analyze_continuous(config, blocks, bad, replicates=0, workers=1)
    bad = deepcopy(parent)
    bad["temperature_K"] = 301
    with pytest.raises(ValueError, match="temperature"):
        module.analyze_continuous(config, blocks, bad, replicates=0, workers=1)


def test_unknown_resampling_modes_and_nonsingleton_iid_groups_fail() -> None:
    options = dict(
        reduced=np.zeros((2, 4)),
        counts=np.array([2, 2]),
        s=np.array([0.2, 0.2, 0.8, 0.8]),
        edges=np.array([0, 0.5, 1]),
        temperature=300,
        reference_s=0.25,
        block_lengths=[[1, 1], [1, 1]],
        replicates=50,
        seed=1,
        workers=1,
    )
    with pytest.raises(ValueError, match="Unknown"):
        block_bootstrap(**options, resampling_unit="unknown")
    options["block_lengths"] = [[2], [2]]
    with pytest.raises(ValueError, match="singleton"):
        block_bootstrap(**options, resampling_unit="decorrelated_configuration")
