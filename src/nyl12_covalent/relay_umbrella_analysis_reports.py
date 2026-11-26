"""Explicit convergence evidence and noninterpolating relay PMF figures."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


def compare_reports(
    current: dict[str, Any], previous: dict[str, Any]
) -> dict[str, Any]:
    """Compare cumulative reports with identical path, binning and E-S reference.

    Args:
        current: Current accepted-block analysis.
        previous: Prior cumulative analysis for exactly the same physical scope.

    Returns:
        Explicit stability evidence, or an unevaluated gate when support is absent.
    """
    for field in (
        "physical_hamiltonian_sha256",
        "temperature_K",
        "projection",
        "path_3cv_sha256",
    ):
        if current[field] != previous[field]:
            raise ValueError(f"Cannot compare reports with different {field}")
    first, second = previous["pmf"], current["pmf"]
    if not first or not second:
        return dict(passed=None, reason="Both reports require connected PMF support")
    for field in ("edges", "es_search_bounds", "acyl_search_bounds"):
        if first[field] != second[field]:
            raise ValueError(f"Cannot compare reports with different {field}")
    if any(item.get("minimum_to_barrier_kcal_mol") is None for item in (first, second)):
        return dict(passed=None, reason="Barrier support is gapped or absent")
    # Path length converts normalized arclength stability into Angstrom.
    length = float(current["path_length_A"])
    barrier_change = abs(
        second["minimum_to_barrier_kcal_mol"] - first["minimum_to_barrier_kcal_mol"]
    )
    location_change = abs(first["barrier_s"] - second["barrier_s"]) * length
    reference_change = abs(first["reference_s"] - second["reference_s"]) * length
    acyl_change = abs(first["acyl_reference_s"] - second["acyl_reference_s"]) * length
    basin_change = abs(
        first["basin_minimum_difference_kcal_mol"]
        - second["basin_minimum_difference_kcal_mol"]
    )
    return dict(
        passed=basin_change <= 0.5
        and barrier_change <= 0.5
        and max(location_change, reference_change, acyl_change) <= 0.1,
        barrier_change_kcal_mol=barrier_change,
        basin_minimum_change_kcal_mol=basin_change,
        barrier_location_change_A=location_change,
        reference_location_change_A=reference_change,
        acyl_location_change_A=acyl_change,
        reason="One cumulative comparison; two successive additions still required",
        promotion_complete=False,
    )


def plot_report(report: dict[str, Any], output: Path) -> None:
    """Save a standalone diagnostic figure without connecting empty PMF bins.

    Args:
        report: Analysis report, including explicit readiness/disconnection status.
        output: New PNG/PDF output path; existing artifacts are never overwritten.
    """
    if output.exists():
        raise FileExistsError(output)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    pmf = report.get("pmf")
    if pmf and pmf.get("pmf_kcal_mol"):
        values = [np.nan if x is None else x for x in pmf["pmf_kcal_mol"]]
        axes[0].plot(pmf["centers"], values, ".-", linewidth=1)
        axes[0].set(
            xlabel="s (raw normalized path arclength)", ylabel="F(s) (kcal/mol)"
        )
        axes[0].set_title("Preliminary conditional modeled channel")
    else:
        axes[0].text(
            0.5,
            0.5,
            "PMF unavailable\n" + report["status"],
            ha="center",
            va="center",
            wrap=True,
        )
        axes[0].set_axis_off()
    if "overlap" in report:
        artist = axes[1].imshow(report["overlap"], origin="lower", vmin=0, vmax=1)
        figure.colorbar(artist, ax=axes[1], label="MBAR overlap")
        axes[1].set(
            xlabel="Window", ylabel="Window", title="Shared chart normalization"
        )
    else:
        axes[1].bar(range(len(report["N_k"])), report["N_k"])
        axes[1].set(
            xlabel="Window", ylabel="Decorrelated N_k", title="Coverage readiness"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def chemical_diagnostics(
    blocks: list[dict[str, Any]], columns: list[str]
) -> dict[str, Any]:
    """Retain owner changes and split-block orthogonal diagnostics for review.

    Args:
        blocks: Accepted records including original uncensored observable rows.
        columns: Numeric CN/orthogonal observable column labels.

    Returns:
        Label-swap flags and split-half summaries. Multimodality and stationarity
        require scientific review; these descriptive metrics do not certify them.
    """
    summaries, changed = [], 0
    owner_states: Counter[str] = Counter()
    for block in blocks:
        rows = block.get("original_rows", [])
        changed += sum(
            str(row.get("owner_changed", "false")).lower() == "true" for row in rows
        )
        owner_states.update(
            "|".join(
                f"{key}={row[key]}"
                for key in sorted(row)
                if key.endswith("_owner") or key == "hydrogen_owners"
            )
            for row in rows
        )
        values = np.asarray(block["values"])
        for index, name in enumerate(columns, start=3):
            first, second = np.array_split(values[:, index], 2)
            summaries.append(
                dict(
                    window=block["window"],
                    observable=name,
                    first_half_mean=float(first.mean()),
                    second_half_mean=float(second.mean()),
                    mean_shift=float(second.mean() - first.mean()),
                    range=[
                        float(values[:, index].min()),
                        float(values[:, index].max()),
                    ],
                )
            )
    return dict(
        owner_changed_frames=changed,
        label_swap_present=changed > 0,
        owner_states=dict(owner_states),
        chemical_review_hold=True,
        identity_or_hidden_multimodality_review="required before promotion; retain every branch",
        orthogonal_stationarity="unresolved",
        split_block_diagnostics=summaries,
        samples_censored=0,
    )
