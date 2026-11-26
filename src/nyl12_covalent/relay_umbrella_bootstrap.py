"""Explicit contiguous-block bootstrap for accepted relay production series."""

from __future__ import annotations

from typing import Any

import numpy as np

from .relay_umbrella_analysis import analyze


def bootstrap(
    windows: list[dict[str, Any]],
    blocks: list[dict[str, Any]],
    *,
    analysis_options: dict[str, Any],
    replicates: int,
    chunk_frames: int,
    seed: int,
) -> dict[str, Any]:
    """Recompute joint MBAR under contiguous moving-block bootstrap resampling.

    Args:
        windows: Complete bias manifest used by the original analysis.
        blocks: Uncensored accepted source-block arrays.
        analysis_options: Exact keyword arguments passed to ``analyze``.
        replicates: At least twenty independently resampled estimates.
        chunk_frames: Explicit correlation-informed contiguous bootstrap length.
            Must exceed the retained series' measured statistical inefficiencies.
        seed: Declared reproducible random seed.

    Returns:
        Percentile 95% intervals and half widths, with failed-support replicates
        retained as failures. Any failed replicate prevents gate promotion.
    """
    if replicates < 20 or chunk_frames < 2:
        raise ValueError(
            "At least twenty replicates and two-frame contiguous chunks required"
        )
    baseline = analyze(windows, blocks, **analysis_options)
    measured = max(
        (item["statistical_inefficiency"] for item in baseline["blocks"]),
        default=np.inf,
    )
    if chunk_frames < measured:
        raise ValueError("Bootstrap chunk length is shorter than measured correlation")
    rng = np.random.default_rng(seed)
    estimates, failures = [], 0
    for _ in range(replicates):
        resampled = []
        for block, summary in zip(blocks, baseline["blocks"], strict=True):
            discard = summary["discarded_frames"]
            values = np.asarray(block["values"])[discard:]
            native = np.asarray(block["native_bias_kcal_mol"])[discard:]
            # Resampling changes frame order; independent native parity was
            # already checked on the immutable baseline, including its terminal mask.
            count = len(values)
            if chunk_frames > count:
                raise ValueError(
                    "Bootstrap chunks cannot exceed the accepted block length"
                )
            starts = rng.integers(
                0, count - chunk_frames + 1, size=int(np.ceil(count / chunk_frames))
            )
            indices = np.concatenate(
                [np.arange(start, start + chunk_frames) for start in starts]
            )[:count]
            resampled.append(
                dict(
                    block,
                    values=values[indices],
                    native_bias_kcal_mol=native[indices],
                    native_bias_evaluated=np.asarray(
                        block.get(
                            "native_bias_evaluated", np.ones(len(block["values"]), bool)
                        )
                    )[discard:][indices],
                )
            )
        for block in resampled:
            block.pop("native_cv_values_A", None)
            block.pop("native_bias_evaluated", None)
            from .relay_umbrella_cv import RelayBiasWindow, cross_bias_energies

            definition = next(
                window["bias"]
                for window in windows
                if window["label"] == block["window"]
            )
            block["native_bias_kcal_mol"] = cross_bias_energies(
                dict(
                    zip(
                        ("q_acyl", "q_act", "q_del"),
                        block["values"][:, :3].T,
                        strict=True,
                    )
                ),
                [RelayBiasWindow(**definition)],
            )[0]
            block["native_bias_source"] = "bootstrap_computed_after_verified_baseline"
        result = analyze(windows, resampled, **analysis_options)
        pmf = result.get("pmf")
        if not pmf or pmf.get("minimum_to_barrier_kcal_mol") is None:
            failures += 1
            continue
        estimates.append(
            [
                pmf["minimum_to_barrier_kcal_mol"],
                pmf["basin_minimum_difference_kcal_mol"],
            ]
        )
    intervals = np.percentile(estimates, [2.5, 97.5], axis=0) if estimates else None
    halfwidth = (intervals[1] - intervals[0]) / 2 if intervals is not None else None
    return dict(
        passed=bool(failures == 0 and halfwidth is not None and np.all(halfwidth <= 1)),
        replicates=replicates,
        failed_support_replicates=failures,
        seed=seed,
        contiguous_chunk_frames=chunk_frames,
        observable_order=[
            "minimum_to_barrier_kcal_mol",
            "basin_minimum_difference_kcal_mol",
        ],
        percentile_95CI=intervals.tolist() if intervals is not None else None,
        halfwidth_kcal_mol=halfwidth.tolist() if halfwidth is not None else None,
        estimates=estimates,
    )
