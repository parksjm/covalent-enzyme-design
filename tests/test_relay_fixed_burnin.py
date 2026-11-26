"""Tests for fixed relay-window burn-in and contiguous chunks."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.relay_fixed_burnin import (
    block_lengths_for,
    concatenate_windows,
    fixed_window_data,
)


def make_block(label: str, start: float, frames: int = 400) -> dict:
    dt = 0.05
    times = start + dt * np.arange(1, frames + 1)
    return dict(
        window=label,
        values=np.column_stack((np.arange(frames), np.arange(frames) ** 2)),
        original_rows=[dict(time_ps=str(time)) for time in times],
        sampled_duration_ps=frames * dt,
    )


def test_join_resets_absolute_start_and_keeps_first_frame_interval() -> None:
    windows = [dict(label="a"), dict(label="b")]
    values, times = concatenate_windows(
        windows, [make_block("a", 100.0), make_block("b", 250.0)]
    )
    assert [len(array) for array in values] == [400, 400]
    assert times[0][0] == pytest.approx(0.05)
    assert times[1][0] == pytest.approx(0.05)
    assert times[0][-1] == pytest.approx(20.0)


def test_fixed_burnin_keeps_contiguous_post_cut_frames_and_summaries() -> None:
    windows = [dict(label="w")]
    samples, counts, summaries = fixed_window_data(
        windows, [make_block("w", 0.0)], burn_in_ps=5.0, end_ps=20.0
    )
    assert counts == [300]
    assert samples[0, 0] == 100
    assert samples[-1, 0] == 399
    summary = summaries[0]
    assert summary["selected_start_ps"] == pytest.approx(5.05)
    assert summary["selected_frames"] == 300
    assert len(summary["g_by_column_diagnostic_only"]) == 2
    assert summary["five_ps_bin_means"][0]["start_ps"] == 5.0
    assert summary["five_ps_bin_means"][-1]["end_ps"] == 20.0
    assert summary["start_bin_mean"] == summary["five_ps_bin_means"][0]["means"]
    assert summary["end_bin_mean"] == summary["five_ps_bin_means"][-1]["means"]


def test_continuity_rejects_duplicate_and_gap_between_segments() -> None:
    windows = [dict(label="w")]
    first = make_block("w", 0.0, frames=100)
    second = make_block("w", 4.95, frames=100)
    with pytest.raises(ValueError, match="Duplicate, gap"):
        concatenate_windows(windows, [first, second])
    second = make_block("w", 5.05, frames=100)
    with pytest.raises(ValueError, match="Duplicate, gap"):
        concatenate_windows(windows, [first, second])


def test_rejects_irregular_spacing_or_unrepresented_end() -> None:
    windows = [dict(label="w")]
    block = make_block("w", 0.0)
    block["original_rows"][10]["time_ps"] = "0.6"
    with pytest.raises(ValueError, match="Irregular frame spacing"):
        concatenate_windows(windows, [block])
    with pytest.raises(ValueError, match="end_ps is not on the frame grid"):
        fixed_window_data(windows, [make_block("w", 0.0)], 5.0, 19.975)


def test_exact_chunk_partition_and_invalid_chunks() -> None:
    assert block_lengths_for(300, 0.05, 2.5) == [50] * 6
    assert block_lengths_for(300, 0.05, 5.0) == [100] * 3
    with pytest.raises(ValueError, match="divide"):
        block_lengths_for(301, 0.05, 2.5)
    with pytest.raises(ValueError, match="integer number"):
        block_lengths_for(300, 0.05, 2.51)


def test_offset_receipt_times_use_integer_grid_for_exact_bounds() -> None:
    windows = [dict(label="offset")]
    samples, counts, summaries = fixed_window_data(
        windows,
        [make_block("offset", 6.0, frames=500)],
        burn_in_ps=5.0,
        end_ps=25.0,
    )
    assert counts.tolist() == [400]
    assert len(samples) == 400
    assert summaries[0]["selected_start_ps"] == pytest.approx(5.05)
    assert summaries[0]["end_ps"] == pytest.approx(25.0)


def test_zero_start_is_a_valid_grid_boundary() -> None:
    samples, counts, _ = fixed_window_data(
        [dict(label="w")],
        [make_block("w", 6.0)],
        burn_in_ps=0.0,
        start_ps=0.0,
        end_ps=20.0,
    )
    assert counts.tolist() == [400]
    assert len(samples) == 400


def test_float32_native_times_preserve_exact_frame_counts() -> None:
    blocks = [make_block("w", 6.0 + 5 * i, frames=100) for i in range(5)]
    for block in blocks:
        for row in block["original_rows"]:
            row["time_ps"] = str(float(np.float32(row["time_ps"])))
    _, counts, summaries = fixed_window_data([dict(label="w")], blocks, 5.0, 25.0)
    assert counts.tolist() == [400]
    assert summaries[0]["frame_interval_ps"] == pytest.approx(0.05)


def test_small_but_cumulative_timestamp_drift_is_rejected() -> None:
    block = make_block("w", 0.0, frames=100)
    for i, row in enumerate(block["original_rows"]):
        row["time_ps"] = str(float(row["time_ps"]) + 0.00002 * i)
    with pytest.raises(ValueError, match="Cumulative timestamp drift"):
        concatenate_windows([dict(label="w")], [block])
