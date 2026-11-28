"""Most-extensive candidate analysis, isolated from all published PMFs."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .relay_final_bootstrap import circular_bootstrap, nonoverlap_baseline
from .relay_final_statistics import summarize
from .relay_pooled_analysis import write
from .relay_pooled_estimator import evaluate, fit_logweights, point, projection
from .relay_umbrella_analysis_io import digest


def definitions(data: dict) -> dict:
    """Build existing BGP sensitivity and explicitly gridded delivery coordinates."""
    from .relay_pooled_estimator import make_projections

    q = data["q"]
    result = make_projections(q, data["config"])
    path = np.asarray(data["config"]["path_3cv_A"])
    nodes, _ = projection.uniform_nodes(path, 153)
    for factor in (0.5, 2.0):
        values, _ = projection.smooth_progress(q, nodes, 56.36835607868901 * factor)
        result[f"BGP_lambda_x{factor}"] = dict(
            values=values,
            edges=np.linspace(0, 1, 41),
            mask=np.ones(len(q), bool),
            reference=0,
        )
    ex, ey = np.arange(-5, 4.00001, 0.1), np.arange(-3, 3.00001, 0.1)
    ix = np.searchsorted(ex, q[:, 2], side="right") - 1
    iy = np.searchsorted(ey, q[:, 0], side="right") - 1
    nx, ny = len(ex) - 1, len(ey) - 1
    if np.any((ix < 0) | (ix >= nx) | (iy < 0) | (iy >= ny)):
        raise ValueError("Frozen delivery grid overflow; no clipping permitted")
    rx = np.searchsorted(ex, path[0, 2], side="right") - 1
    ry = np.searchsorted(ey, path[0, 0], side="right") - 1
    result["delivery_2D"] = dict(
        values=ix * ny + iy + 0.5,
        edges=np.arange(nx * ny + 1),
        mask=np.ones(len(q), bool),
        reference=int(rx * ny + ry),
        x_edges=ex,
        y_edges=ey,
        shape=(nx, ny),
    )
    return result


def plot(report: dict, grids: dict, output: Path) -> None:
    """Plot the actual final candidate without any primary alias or smoothing."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for name in ("hard", "BGP", "BGP_lambda_x0.5", "BGP_lambda_x2.0"):
        p = grids[name]
        axes[0, 0].plot(
            (p["edges"][:-1] + p["edges"][1:]) / 2,
            report["projections"][name]["pmf"],
            label=name,
        )
    axes[0, 0].legend()
    axes[0, 0].set(
        xlabel="Named progress coordinate", ylabel="Own-reference ΔF [kcal/mol]"
    )
    for ax, name, xlabel in [
        (axes[0, 1], "initial_2D", "q_act [Å]"),
        (axes[1, 0], "delivery_2D", "q_del [Å]"),
    ]:
        p = grids[name]
        f = np.asarray(report["projections"][name]["pmf"], float).reshape(p["shape"])
        mesh = ax.pcolormesh(p["x_edges"], p["y_edges"], f.T, vmin=-2, vmax=30)
        ax.set(xlabel=xlabel, ylabel="q_acyl [Å]", title=name + "; global MBAR weights")
        fig.colorbar(mesh, ax=ax, label="Reference-cell ΔF [kcal/mol]")
    h = grids["hard"]
    axes[1, 1].semilogy(
        (h["edges"][:-1] + h["edges"][1:]) / 2, report["projections"]["hard"]["ESS"]
    )
    axes[1, 1].set(xlabel="Hard s", ylabel="Kish weight ESS; not temporal ESS")
    fig.suptitle(
        "FINAL CANDIDATE: main57 actual60ps +43 independent30ps lineages\n"
        "Not promoted; stationarity unresolved; no cosmetic KDE"
    )
    fig.savefig(output / "final_candidate_point.png", dpi=170)
    plt.close(fig)


def temporal_halves(reduced: np.ndarray, counts: np.ndarray, grids: dict) -> dict:
    """Refit matched first/last retained halves independently in every lineage."""
    result = {}
    offsets = np.r_[0, np.cumsum(counts)]
    for half in (0, 1):
        ids = np.concatenate(
            [
                np.arange(
                    a + (0 if half == 0 else n // 2), a + (n // 2 if half == 0 else n)
                )
                for a, n in zip(offsets[:-1], counts, strict=True)
            ]
        )
        half_counts = counts // 2 if half == 0 else counts - counts // 2
        lw, fit = fit_logweights(reduced[:, ids], half_counts)
        result["first" if half == 0 else "last"] = dict(
            connected=fit["connected"],
            counts=half_counts.tolist(),
            projections=evaluate(lw, grids, ids),
            interpretation="Matched retained halves; temporal drift diagnostic, not independent replicas",
        )
    return result


def analyze(data: dict, output: Path, workers: int = 4) -> dict:
    """Run one global point fit and sequential paired 5/10ps uncertainty methods."""
    output.mkdir(parents=True, exist_ok=False)
    report, _, reduced, lw = point(data)
    grids = definitions(data)
    report["projections"] = evaluate(lw, grids, np.arange(len(lw)))
    report.update(
        main_actual_horizon_ps=60,
        supplemental_horizon_ps=30,
        CI_gate="Pending 200 circular5/10ps and nonoverlap5ps refits",
        warning="A026–028 nonstationary; gdt up to11.8ps. Block CIs do not establish mixing.",
        sampling_converged=False,
        primary_replacement=False,
    )
    write(output / "point.json", report)
    write(
        output / "sources.json",
        {k: v for k, v in data.items() if k not in ("q", "counts")},
    )
    np.savez_compressed(
        output / "matched_final_weights.npz",
        q3cv_A=data["q"],
        counts=data["counts"],
        log_weights=lw,
        reduced_bias=reduced,
    )
    plot(report, grids, output)
    write(output / "physical_lineage_statistics.json", summarize(data))
    write(
        output / "temporal_half_fits.json",
        temporal_halves(reduced, data["counts"], grids),
    )
    for block in (100, 200):
        result = circular_bootstrap(reduced, data["counts"], grids, block, workers)
        write(output / f"circular_{block // 20}ps_conditional95.json", result)
    write(
        output / "nonoverlap_5ps_comparator.json",
        nonoverlap_baseline(reduced, data["counts"], grids, workers),
    )
    acceptance = dict(
        software_complete=True,
        states=100,
        frames=len(lw),
        actual_main_ps=60,
        actual_supplement_ps=30,
        sampling_converged=False,
        primary_replacement=False,
        point_sha256=digest(output / "point.json"),
    )
    write(output / "acceptance.json", acceptance)
    return acceptance
