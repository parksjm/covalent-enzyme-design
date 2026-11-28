"""Isolated pooled-support reports and conditional uncertainty; never primary aliases."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .relay_pooled_estimator import nullable, paired_bootstrap, point
from .relay_umbrella_analysis_io import digest


def serializable(value: object) -> object:
    """Convert numerical artifacts to authoritative JSON-compatible values."""
    if isinstance(value, np.ndarray):
        return nullable(value)
    if isinstance(value, dict):
        return {k: serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write(path: Path, value: dict) -> None:
    """Retain a complete JSON diagnostic without NaN-to-zero substitution."""
    path.write_text(json.dumps(serializable(value), indent=2, allow_nan=False) + "\n")


def plot_report(
    report: dict, definitions: dict, output: Path, intervals: dict | None = None
) -> None:
    """Write legible hard/BGP/physical/2D diagnostics without altering primary plots."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for name, color in [("hard", "black"), ("BGP", "tab:blue")]:
        p = definitions[name]
        h = report["projections"][name]
        centers = (p["edges"][:-1] + p["edges"][1:]) / 2
        axes[0, 0].plot(centers, h["pmf"], label=name, color=color)
        if intervals is not None:
            ci = intervals[name]
            axes[0, 0].fill_between(
                centers,
                np.asarray(ci["lower_95"], float),
                np.asarray(ci["upper_95"], float),
                alpha=0.18,
                color=color,
            )
    axes[0, 0].set(
        xlabel="Named progress coordinate (different measures)",
        ylabel="Own-reference ΔF [kcal/mol]",
        title="Pooled global point fit; original primary unchanged",
    )
    axes[0, 0].legend()
    name = "initial_q_act"
    p = definitions[name]
    h = report["projections"][name]
    centers = (p["edges"][:-1] + p["edges"][1:]) / 2
    axes[0, 1].plot(centers, h["pmf"])
    if intervals is not None:
        axes[0, 1].fill_between(
            centers,
            np.asarray(intervals[name]["lower_95"], float),
            np.asarray(intervals[name]["upper_95"], float),
            alpha=0.2,
        )
    axes[0, 1].set(
        xlim=(-1.4, 1.4),
        xlabel="Physical q_act [Å]",
        ylabel="Conditional ΔF [kcal/mol]",
        title="Initial branch −2≤q_acyl<−0.8Å",
    )
    p = definitions["initial_2D"]
    f = np.array(report["projections"]["initial_2D"]["pmf"], float).reshape(p["shape"])
    mesh = axes[1, 0].pcolormesh(
        p["x_edges"],
        p["y_edges"],
        f.T,
        vmin=-2,
        vmax=20,
        shading="flat",
        cmap="viridis",
    )
    axes[1, 0].set(
        xlim=(-1.4, 1.4),
        xlabel="q_act [Å]",
        ylabel="q_acyl [Å]",
        title="Initial-branch 2D; same global weights",
    )
    fig.colorbar(mesh, ax=axes[1, 0], label="Reference-cell ΔF [kcal/mol]")
    hard = definitions["hard"]
    centers = (hard["edges"][:-1] + hard["edges"][1:]) / 2
    axes[1, 1].semilogy(
        centers, report["projections"]["hard"]["ESS"], label="Weight ESS"
    )
    axes[1, 1].set(
        xlabel="Hard progress s",
        ylabel="Kish weight ESS (not temporal ESS)",
        title="Support, including product tail",
    )
    subtitle = (
        "POINT ONLY; CIs UNRESOLVED/pending; see strict source chunk gate"
        if intervals is None
        else "Pointwise conditional95%; same200lineage-block resamples for all panels"
    )
    fig.suptitle(
        f"{report['states']} independent MBAR states · main57 frozen45ps · supplement{report['supplemental_horizon_ps']}ps\n{subtitle}"
    )
    fig.savefig(
        output
        / ("pooled_point.png" if intervals is None else "pooled_conditional_CI.png"),
        dpi=170,
    )
    plt.close(fig)


def analyze(data: dict, output: Path, workers: int = 4) -> dict:
    """Fit one isolated checkpoint and eligible paired block uncertainty.

    Args:
        data: Explicit source-bound independent-lineage observations.
        output: New diagnostic checkpoint directory.
        workers: At most4bootstrapworkers.

    Returns:
        Point/uncertainty acceptance record; publication always requires lead review.
    """
    if not 1 <= workers <= 4:
        raise ValueError("CPUlimitfourworkers")
    output.mkdir(parents=True, exist_ok=False)
    source = {k: v for k, v in data.items() if k not in ("q", "counts")}
    write(output / "sources.json", source)
    report, definitions, reduced, lw = point(data)
    report["supplemental_horizon_ps"] = data["supplemental_horizon_ps"]
    report["comparison_scope"] = (
        "CONTROLLED frozen-main45ps supplemental-support diagnostic; NOT final or primary replacement"
    )
    report["final_most_extensive_analysis_required"] = (
        "All accepted main60ps + oldtarget30ps + new33_30ps; independent5psdiscard and pairedbootstrap"
    )
    report["lineages"] = data["states"]
    report["stationarity_warning"] = (
        "Activation026-028 show HG1-Nleave mean drift ~1.53-1.62A and zero early/late histogram overlap; C-OG1 also drifts. 5/2.5ps bootstrap CIs can understate nonstationarity and are NOT convergence validation. Early33grid does not by itself fix this region."
    )
    report["final_validation_required"] = (
        "Most extensive main60ps + old4_30ps + new33_30ps, with >=10ps block sensitivity or independent-start evidence; never assume five-ps blocks independent."
    )
    np.savez_compressed(
        output / "global_weights.npz",
        q=data["q"],
        counts=data["counts"],
        logweights=lw,
        f_k=report["f_k"],
        state_index=np.repeat(np.arange(len(data["counts"])), data["counts"]),
    )
    write(output / "point.json", report)
    plot_report(report, definitions, output)
    write(
        output / "projection_definitions.json",
        {
            k: {x: v for x, v in p.items() if x not in ("values", "mask")}
            for k, p in definitions.items()
        },
    )
    if not data["bootstrap_eligible"]:
        write(
            output / "acceptance.json",
            dict(
                point_fit=True,
                states=report["states"],
                horizon_ps=data["supplemental_horizon_ps"],
                CI_status="UNRESOLVED:onlyone5pschunkpernewindependentlineage;no2.5psfallback",
                bootstrap_refits=0,
                primary_replacement=False,
                source_sha256=digest(output / "sources.json"),
            ),
        )
        return report
    intervals = {}
    for width, frames in [("5ps", 100), ("2.5ps", 50)]:
        intervals[width] = paired_bootstrap(
            reduced, data["counts"], definitions, frames, workers=workers
        )
        write(output / f"bootstrap_{width}.json", intervals[width])
    # Matched within-lineage halves assess temporal drift separately from CIs.
    half_reports = {}
    for name, late in [("early", False), ("late", True)]:
        selected = []
        offset = 0
        for count in data["counts"]:
            half = int(count) // 2
            selected.extend(
                range(
                    offset + half if late else offset,
                    offset + int(count) if late else offset + half,
                )
            )
            offset += int(count)
        ids = np.asarray(selected)
        half = dict(data, q=data["q"][ids], counts=data["counts"] // 2)
        half_reports[name] = point(half)[0]
    write(output / "temporal_halves.json", half_reports)
    plot_report(report, definitions, output, intervals["5ps"])
    write(
        output / "acceptance.json",
        dict(
            point_fit=True,
            states=report["states"],
            horizon_ps=data["supplemental_horizon_ps"],
            bootstrap_refits=400,
            replicates_per_width=200,
            workers=workers,
            CI_status="Pointwiseconditional95%;supportmustbe1forreportedinterval",
            primary_replacement=False,
            sampling_converged=False,
            source_sha256=digest(output / "sources.json"),
        ),
    )
    return report
