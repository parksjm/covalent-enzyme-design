"""Fixed per-window burn-in selection for continuous relay production."""

from __future__ import annotations

from typing import Any

import numpy as np
from pymbar import timeseries


def concatenate_windows(
    windows: list[dict[str, Any]], blocks: list[dict[str, Any]]
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Join accepted production segments and establish a common time grid.

    Each window is ordered by its original timestamps. Absolute start times may
    differ by window; returned times start at one frame interval and increase
    from zero independently for each window.

    Args:
        windows: Ordered umbrella-window definitions containing unique labels.
        blocks: Accepted blocks with ``window``, ``values``, ``original_rows``
            and ``sampled_duration_ps`` fields.

    Returns:
        Window-ordered values and zero-relative production times in picoseconds.

    Raises:
        ValueError: If segments are absent, malformed, discontinuous, or use
            different frame intervals.
    """
    labels = [window["label"] for window in windows]
    if len(labels) < 1 or len(set(labels)) != len(labels):
        raise ValueError("At least one window with unique labels is required")
    if any(block["window"] not in labels for block in blocks):
        raise ValueError("Production block references an unknown window")

    joined_values: list[np.ndarray] = []
    joined_times: list[np.ndarray] = []
    common_dt: float | None = None
    for label in labels:
        parts = [block for block in blocks if block["window"] == label]
        if not parts:
            raise ValueError(f"No production blocks for {label}")
        parts.sort(key=lambda block: float(block["original_rows"][0]["time_ps"]))
        times_parts: list[np.ndarray] = []
        values_parts: list[np.ndarray] = []
        for part in parts:
            values = np.asarray(part["values"], dtype=float)
            rows = part["original_rows"]
            times = np.asarray([float(row["time_ps"]) for row in rows])
            if (
                values.ndim != 2
                or len(values) != len(times)
                or len(times) < 2
                or not np.isfinite(values).all()
                or not np.isfinite(times).all()
                or np.any(np.diff(times) <= 0)
            ):
                raise ValueError(f"Invalid time or observable data for {label}")
            duration = float(part["sampled_duration_ps"])
            if not np.isfinite(duration) or duration <= 0:
                raise ValueError(f"Invalid sampled duration in {label}")
            dt = duration / len(times)
            # Accepted NetCDF times are float32. Recover the exact spacing from
            # duration/count, while checking serialized times to0.0001ps.
            if not np.allclose(np.diff(times), dt, rtol=0, atol=1e-4):
                raise ValueError(f"Irregular frame spacing in {label}")
            if common_dt is None:
                common_dt = dt
            elif not np.isclose(dt, common_dt, rtol=0, atol=1e-6):
                raise ValueError("All windows must use a common frame interval")
            times_parts.append(times)
            values_parts.append(values)
        times = np.concatenate(times_parts)
        values = np.concatenate(values_parts)
        spacing = np.diff(times)
        if not np.allclose(spacing, common_dt, rtol=0, atol=1e-4):
            raise ValueError(f"Duplicate, gap, or discontinuity in {label}")
        expected = times[0] + np.arange(len(times)) * common_dt
        if not np.allclose(times, expected, rtol=0, atol=1e-4):
            raise ValueError(f"Cumulative timestamp drift in {label}")
        relative_times = np.arange(1, len(times) + 1, dtype=float) * common_dt
        joined_values.append(values)
        joined_times.append(relative_times)
    return joined_values, joined_times


def fixed_window_data(
    windows: list[dict[str, Any]],
    blocks: list[dict[str, Any]],
    burn_in_ps: float,
    end_ps: float,
    start_ps: float | None = None,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Select a fixed production interval independently in every window.

    Frames at or before ``burn_in_ps`` are discarded. The selected interval is
    otherwise exact: frames must lie strictly after the optional start and at
    or before ``end_ps``; bounds are never silently clipped.

    Args:
        windows: Ordered umbrella-window definitions containing labels.
        blocks: Accepted production blocks understood by
            :func:`concatenate_windows`.
        burn_in_ps: Fixed per-window burn-in boundary in picoseconds.
        end_ps: Inclusive production end boundary in picoseconds.
        start_ps: Optional strict lower bound for a sensitivity horizon.

    Returns:
        Selected value arrays, per-window sample counts, and JSON-ready
        diagnostics including correlation estimates and 5-ps bin means.

    Raises:
        ValueError: If bounds are invalid or the exact requested interval lacks
            frames in any window.
    """
    bounds = [burn_in_ps, end_ps] + ([] if start_ps is None else [start_ps])
    if not np.isfinite(bounds).all() or burn_in_ps < 0 or end_ps <= burn_in_ps:
        raise ValueError("Finite ordered nonnegative burn-in/end bounds required")
    if start_ps is not None and not burn_in_ps <= start_ps < end_ps:
        raise ValueError("start_ps must be at or after burn-in and before end")

    values_by_window, times_by_window = concatenate_windows(windows, blocks)
    selected: list[np.ndarray] = []
    counts: list[int] = []
    summaries: list[dict[str, Any]] = []
    for window, values, times in zip(
        windows, values_by_window, times_by_window, strict=True
    ):
        lower = burn_in_ps if start_ps is None else max(burn_in_ps, start_ps)
        dt = float(times[1] - times[0])
        lower_index = _frame_index(lower, dt, "lower selection bound")
        end_index = _frame_index(end_ps, dt, "end_ps")
        if end_index > len(times):
            raise ValueError(
                f"Production does not reach exact end_ps for {window['label']}"
            )
        frame_indices = np.arange(1, len(times) + 1)
        mask = (frame_indices > lower_index) & (frame_indices <= end_index)
        chosen = values[mask]
        chosen_times = times[mask]
        if not len(chosen):
            raise ValueError(
                f"No frames within requested interval for {window['label']}"
            )

        inefficiencies = []
        for column in chosen.T:
            inefficiencies.append(
                float(timeseries.statistical_inefficiency(column))
                if np.ptp(column) > 1e-12 and len(column) >= 4
                else 1.0
            )
        bin_means = []
        left = float(lower)
        while left < end_ps:
            right = min(left + 5.0, end_ps)
            bin_mask = (chosen_times > left) & (chosen_times <= right)
            if bin_mask.any():
                bin_means.append(
                    dict(
                        start_ps=left,
                        end_ps=right,
                        frame_count=int(bin_mask.sum()),
                        means=chosen[bin_mask].mean(axis=0).tolist(),
                    )
                )
            left = right
        summary = dict(
            window=window["label"],
            burn_in_ps=burn_in_ps,
            selected_start_ps=float(chosen_times[0]),
            end_ps=float(chosen_times[-1]),
            frame_interval_ps=dt,
            selected_frames=len(chosen),
            retained_frames=len(chosen),
            g_by_column=inefficiencies,
            selected_indices=np.flatnonzero(mask).tolist(),
            g_by_column_diagnostic_only=inefficiencies,
            five_ps_bin_means=bin_means,
            start_bin_mean=bin_means[0]["means"] if bin_means else None,
            end_bin_mean=bin_means[-1]["means"] if bin_means else None,
        )
        selected.append(chosen)
        counts.append(len(chosen))
        summaries.append(summary)
    return np.concatenate(selected), np.asarray(counts, dtype=int), summaries


def block_lengths_for(count: int, dt: float, chunk_ps: float) -> list[int]:
    """Return exact equal frame partitions for a requested block duration.

    Args:
        count: Number of selected contiguous frames.
        dt: Common frame interval in picoseconds.
        chunk_ps: Desired contiguous bootstrap chunk duration in picoseconds.

    Returns:
        Frame counts for equal contiguous chunks.

    Raises:
        ValueError: If units do not yield an exact partition with at least two
            complete chunks.
    """
    if count <= 0 or not np.isfinite([dt, chunk_ps]).all() or dt <= 0 or chunk_ps <= 0:
        raise ValueError("Positive count, frame interval, and chunk duration required")
    frames = chunk_ps / dt
    if not np.isclose(frames, round(frames), rtol=0, atol=1e-8):
        raise ValueError("Chunk duration must be an integer number of frames")
    frames_per_chunk = int(round(frames))
    if count < 2 * frames_per_chunk or count % frames_per_chunk:
        raise ValueError("Selected frames must divide into at least two full chunks")
    return [frames_per_chunk] * (count // frames_per_chunk)


def _frame_index(bound_ps: float, dt: float, name: str) -> int:
    """Validate a time boundary and return its one-based frame index."""
    frame_position = bound_ps / dt
    nearest = int(round(frame_position))
    if not np.isclose(frame_position, nearest, rtol=0, atol=1e-8):
        raise ValueError(f"{name} is not on the frame grid")
    return nearest
