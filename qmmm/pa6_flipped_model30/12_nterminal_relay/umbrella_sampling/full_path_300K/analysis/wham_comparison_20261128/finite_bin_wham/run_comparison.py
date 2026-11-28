"""Run one frozen-data finite-bin WHAM sensitivity analysis, without new sampling."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import logging
from pathlib import Path
import sys
import time

import numpy as np
from scipy.special import logsumexp

from wham_core import projected_pmf, quantize, solve_wham
from wham_figures import plot_comparison
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from nyl12_covalent.relay_umbrella_projection import project_path

KT = 0.00198720425864083 * 300
LOGGER = logging.getLogger(__name__)


def native_bias(q: np.ndarray, windows: list) -> np.ndarray:
    """Return all exact native reduced biases at the supplied physical CVs."""
    return (
        cross_bias_energies(
            dict(zip(("q_acyl", "q_act", "q_del"), q.T, strict=True)), windows
        )
        / KT
    )


def nullable(values: np.ndarray) -> list:
    """Encode unsupported bins as JSON null, preserving numerical signs."""
    return [float(v) if np.isfinite(v) else None for v in values]


def statistics(pmf: np.ndarray, baseline: np.ndarray, s: np.ndarray) -> dict:
    """Compare on the fixed original supported bins in the declared early range."""
    common = np.isfinite(baseline) & (s >= 0) & (s <= 0.3)
    if not common.any() or not np.isfinite(pmf[common]).all():
        raise ValueError("Missing frozen early common support")
    delta = pmf[common] - baseline[common]
    peak = int(np.nanargmax(pmf))
    return dict(
        maximum_occupied_bin_kcal_mol=float(pmf[peak]),
        maximum_occupied_bin_s=float(s[peak]),
        maximum_occupied_bin_index=peak,
        early_bins=np.flatnonzero(common).tolist(),
        early_RMSD_kcal_mol=float(np.sqrt(np.mean(delta**2))),
        early_max_abs_deviation_kcal_mol=float(abs(delta).max()),
    )


def main() -> None:
    """Run frozen-input diagnostic and record independent scientific acceptance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    args.output.mkdir(parents=True, exist_ok=False)
    frozen = np.load(args.input / "weights.npz")
    sources = json.loads((args.input / "sources.json").read_text())
    point = json.loads((args.input / "point.json").read_text())
    q, counts, baseline_lw, native_f = (
        frozen[k] for k in ("q", "counts", "logweights", "f_k")
    )
    if q.shape != (82000, 3) or counts.shape != (100,) or counts.sum() != len(q):
        raise ValueError(
            "Diagnostic is pinned to the approved 82000-frame/100-history pool"
        )
    if sources["config"]["temperature_K"] != 300:
        raise ValueError("Unexpected physical temperature")
    windows = [RelayBiasWindow(**state["bias"]) for state in sources["states"]]
    edges = np.array(sources["config"]["s_bin_edges"])
    path = np.array(sources["config"]["path_3cv_A"])
    s = project_path(q, path)["s"]
    bin_centers = (edges[1:] + edges[:-1]) / 2
    reference = int(np.argmin(abs(bin_centers - 0.0125)))
    baseline, baseline_mass = projected_pmf(s, baseline_lw, edges, reference, KT)
    point_pmf = np.array(point["projections"]["hard"]["pmf"], float)
    if not np.array_equal(np.isfinite(baseline), np.isfinite(point_pmf)):
        raise ValueError("Frozen MBAR support parity failure")
    parity = float(np.nanmax(abs(baseline - point_pmf)))
    if parity > 1e-8:
        raise ValueError(f"Frozen MBAR point parity failure {parity}")
    exact_u = native_bias(q, windows)
    baseline_denominator = logsumexp(
        np.log(counts)[:, None] + native_f[:, None] - exact_u, axis=0
    )
    reconstructed = -baseline_denominator
    reconstructed -= logsumexp(reconstructed)
    weight_parity = float(abs(reconstructed - baseline_lw).max())
    if weight_parity > 1e-7:
        raise ValueError(f"Frozen MBAR weight parity failure {weight_parity}")
    owner = np.repeat(np.arange(len(counts)), counts)
    own_exact = exact_u[owner, np.arange(len(q))]
    del exact_u
    rows = [
        dict(
            label="MBAR (original frame biases)",
            **statistics(baseline, baseline, bin_centers),
        )
    ]
    curves = {"MBAR (original frame biases)": baseline}
    center_curves = {}
    records = []
    for width, shift in [(0.1, 0), (0.05, 0), (0.025, 0), (0.0125, 0), (0.0125, 0.5)]:
        start = time.monotonic()
        label = f"WHAM {width:g} A" + ("; half-cell shift" if shift else "")
        LOGGER.info("Starting %s", label)
        centers, inverse, occupancy = quantize(q, width, shift)
        cell_u = native_bias(centers, windows)
        solved = solve_wham(cell_u, inverse, counts, native_f)
        pmf, mass = projected_pmf(s, solved["logweights"], edges, reference, KT)
        center_s = project_path(centers, path)["s"]
        center_pmf, center_mass = projected_pmf(
            center_s, solved["cell_logmass"], edges, reference, KT
        )
        own_error = abs(cell_u[owner, inverse] - own_exact)
        if not np.isclose(mass.sum(), 1, atol=1e-12, rtol=0):
            raise RuntimeError("Projected WHAM mass normalization failed")
        cell_mass = np.bincount(inverse, weights=np.exp(solved["logweights"]))
        mass_error = float(abs(cell_mass - np.exp(solved["cell_logmass"])).max())
        if mass_error > 1e-12:
            raise RuntimeError("Cell/frame mass mismatch")
        record = dict(
            label=label,
            width_A=width,
            shift_fraction_each_axis=shift,
            occupied_cells=len(centers),
            fixed_point_residual_kBT=solved["fixed_point_residual_kBT"],
            sampled_own_bias_absolute_error_p95_kBT=float(np.percentile(own_error, 95)),
            sampled_own_bias_absolute_error_max_kBT=float(own_error.max()),
            cell_frame_mass_max_error=mass_error,
            elapsed_seconds=time.monotonic() - start,
            **statistics(pmf, baseline, bin_centers),
            original_frame_projection_pmf_kcal_mol=nullable(pmf),
            original_frame_projection_mass=mass.tolist(),
            cell_center_projection_pmf_kcal_mol=nullable(center_pmf),
            cell_center_projection_mass=center_mass.tolist(),
        )
        records.append(record)
        rows.append({k: record[k] for k in rows[0]})
        curves[label], center_curves[label] = pmf, center_pmf
        name = f"width_{width:g}_shift_{shift:g}"
        np.savez_compressed(
            args.output / f"{name}.npz",
            centers_A=centers,
            inverse=inverse,
            occupancy=occupancy,
            logweights=solved["logweights"],
            cell_logmass=solved["cell_logmass"],
            f_k=solved["f_k"],
        )
        (args.output / f"{name}.json").write_text(json.dumps(record, indent=2) + "\n")
        LOGGER.info(
            "Completed %s: %s",
            label,
            {k: v for k, v in record.items() if not isinstance(v, list)},
        )
    finest = curves["WHAM 0.0125 A"]
    previous = curves["WHAM 0.025 A"]
    supported = np.isfinite(baseline) & (bin_centers >= 0) & (bin_centers <= 0.3)
    shift_curve = curves["WHAM 0.0125 A; half-cell shift"]
    result = dict(
        estimator="Finite-bin WHAM solved by PyMBAR; shared solver, no independent-software claim",
        projection="Original-frame hard s, equal WHAM frame weights within physical-q cells",
        center_projection="Separate discrete center marginal; centers can cross s-bin boundaries",
        widths_A=[0.1, 0.05, 0.025, 0.0125],
        anchor_A=[0, 0, 0],
        shifted_grid_fraction=[0.5, 0.5, 0.5],
        reference_bin=reference,
        reference_s=float(bin_centers[reference]),
        s_edges=edges.tolist(),
        baseline_PMF_kcal_mol=nullable(baseline),
        baseline_mass=baseline_mass.tolist(),
        baseline_statistics=rows[0],
        frames=len(q),
        histories=len(counts),
        unique_biases=len(
            {json.dumps(st["bias"], sort_keys=True) for st in sources["states"]}
        ),
        retained_counts=counts.tolist(),
        baseline_point_parity_max_kcal_mol=parity,
        baseline_logweight_parity_max=weight_parity,
        finest_vs_previous_early_max_kcal_mol=float(
            abs(finest[supported] - previous[supported]).max()
        ),
        finest_halfshift_early_max_kcal_mol=float(
            abs(finest[supported] - shift_curve[supported]).max()
        ),
        records=records,
        uncertainty="No new confidence intervals",
        sampling_converged=False,
        independent_solver_crosscheck=False,
        versions={
            k: importlib.metadata.version(k)
            for k in ("numpy", "scipy", "pymbar", "matplotlib")
        },
        python=sys.executable,
        input_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(args.input.iterdir())
        },
    )
    (args.output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    with (args.output / "comparison.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot_comparison(bin_centers, curves, center_curves, args.output)
    acceptance = dict(
        process_completion="See separate launcher status.json",
        scientific_checks_passed=True,
        checks=[
            "Frozen MBAR PMF and weight parity",
            "All 100 histories and 82000 frames retained",
            "Native 3-CV bias evaluated without dropping restraints",
            "WHAM fixed-point residual <=1e-7 kBT",
            "Cell/frame mass agreement <=1e-12",
            "Original observable bins/reference preserved",
        ],
        interpretation="Estimator discretization sensitivity only; no sampling/chemical convergence or activation barrier acceptance",
        sampling_converged=False,
        new_simulation=False,
    )
    (args.output / "scientific_acceptance.json").write_text(
        json.dumps(acceptance, indent=2) + "\n"
    )
    LOGGER.info("All five finite-bin WHAM comparisons and numerical checks completed")


if __name__ == "__main__":
    main()
