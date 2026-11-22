"""Reproduce native bias MBAR and an explicit projected conditional path PMF."""
from __future__ import annotations

import csv
import importlib.metadata
import json
import logging
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pymbar import MBAR
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_neb_umbrella import refresh
from nyl12_covalent.neb_umbrella_analysis import KB, collect, components
from nyl12_covalent.neb_umbrella_bias import HarmonicWindow, bias_energy


def project(points: np.ndarray, centers: np.ndarray) -> np.ndarray:
    """Project each native two-CV point onto nearest polyline segment arc length."""
    delta = np.diff(centers, axis=0)
    lengths = np.linalg.norm(delta, axis=1)
    if np.any(lengths <= 0):
        raise ValueError("Duplicate consecutive path centers")
    t = np.clip(np.einsum("nsi,si->ns", points[:, None] - centers[:-1], delta)
                / lengths**2, 0, 1)
    residual = points[:, None] - (centers[:-1] + t[:, :, None] * delta)
    nearest = np.argmin(np.sum(residual**2, axis=2), axis=1)
    return (np.r_[0, np.cumsum(lengths)][nearest]
            + t[np.arange(len(points)), nearest] * lengths[nearest])


def solve(samples: list[np.ndarray], biases: list[HarmonicWindow], beta: float,
          initial: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, bool]:
    """Fit MBAR with the campaign's exact piecewise native Amber bias function."""
    points = np.concatenate(samples)
    reduced = np.array([[beta * bias_energy(tuple(p[:2]), b) for p in points]
                        for b in biases])
    counts = np.array([len(s) for s in samples])
    mbar = MBAR(reduced, counts, verbose=False, solver_protocol="robust",
                initial_f_k=initial)
    logw = -logsumexp(np.log(counts)[:, None] + mbar.f_k[:, None] - reduced, axis=0)
    weights = np.exp(logw - logsumexp(logw))
    connected = len(components(mbar.compute_overlap()["matrix"])) == 1
    return weights, mbar.f_k, connected


def profile(s: np.ndarray, weights: np.ndarray, edges: np.ndarray,
            beta: float) -> np.ndarray:
    """Return equal-width path marginal relative to the fixed first ES bin."""
    if np.any(s < edges[0]) or np.any(s > edges[-1]):
        raise ValueError("Projected coordinate outside declared path bounds")
    indices = np.minimum(np.searchsorted(edges, s, side="right") - 1, len(edges) - 2)
    hist = np.bincount(indices, weights=weights, minlength=len(edges) - 1)
    answer = np.full(len(hist), np.nan)
    valid = hist > 0
    answer[valid] = -np.log(hist[valid] / np.diff(edges)[valid]) / beta
    if not valid[0]:
        raise ValueError("ES reference bin empty")
    return answer - answer[0]


def main() -> None:
    """Run frozen MBAR, block-bootstrap path marginal, plots and acceptance record."""
    logging.basicConfig(level=logging.INFO)
    out = ROOT.parent
    cfg, windows = collect(out / "frozen_inputs/prepared.json")
    result = refresh(out / "frozen_inputs/prepared.json", out / "mbar",
                     plot_python=Path(sys.executable))
    if not result["connected"]:
        raise ValueError("Disconnected sampled overlap; no global projected PMF")
    centers = np.array([[w["q_acyl_center"], w["q_h_center"]]
                        for w in cfg["windows"]])
    biases = [HarmonicWindow(**w["bias"]) for w in cfg["windows"]]
    strides = [w["decorrelation_stride"] for w in result["window_summaries"]]
    raw = [np.concatenate([b["samples"] for b in w["blocks"]]) for w in windows]
    selected = [s[::stride] for s, stride in zip(raw, strides, strict=True)]
    beta = 1 / (KB * cfg["temperature_K"])
    weights, free, connected = solve(selected, biases, beta)
    assert connected
    assert np.isclose(1 / np.sum(weights**2), result["unbiased_effective_sample_size"])
    edges = np.linspace(0, np.linalg.norm(np.diff(centers, axis=0), axis=1).sum(), 26)
    samples = np.concatenate(selected)
    coordinate = project(samples[:, :2], centers)
    pmf = profile(coordinate, weights, edges, beta)
    rng = np.random.default_rng(20261121)
    bootstrap = []
    disconnected = 0
    for replicate in range(50):
        resampled = []
        for w, stride in zip(windows, strides, strict=True):
            blocks = w["blocks"]
            picks = rng.integers(0, len(blocks), size=len(blocks))
            resampled.append(np.concatenate([blocks[i]["samples"] for i in picks])[::stride])
        bw, _, ok = solve(resampled, biases, beta, free)
        if not ok:
            disconnected += 1
            continue
        bs = project(np.concatenate(resampled)[:, :2], centers)
        bootstrap.append(profile(bs, bw, edges, beta))
        logging.info("Bootstrap %d/50", replicate + 1)
    boot = np.array(bootstrap)
    low, high = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    mid = (edges[:-1] + edges[1:]) / 2
    maximum = int(np.nanargmax(pmf))
    with (out / "path_pmf.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["s_arc_A", "F_relative_ES_kcal_mol", "block_bootstrap_low95", "block_bootstrap_high95"])
        writer.writerows(zip(mid, pmf, low, high, strict=True))
    block_counts = [len(w["blocks"]) for w in windows]
    summaries = result["window_summaries"]
    summary = dict(
        scientific_identity="300K DFTB3/3ob, historical neutral-Thr-N ASP272/ASP274, fixed exterior beyond12A; conditional preliminary free energy",
        projection="Nearest orthogonal segment projection onto ordered (q_acyl,q_H) umbrella-center polyline; Euclidean arc length in Angstrom; transverse coordinates integrated, no tube cutoff",
        reference="First equal-width path bin adjoining ES; bin density reference, not an integrated ES basin",
        exact_bias="Imported unchanged native Amber k*d^2 plus tangent linear tails; both chemical CVs included in every cross-bias",
        uncertainty="Within-window resampling of completed5ps blocks with replacement, MBAR refit each replicate, original per-window statistical-thinning stride fixed; percentile95 interval; finite-sample exploratory only, not convergence proof",
        bootstrap_seed=20261121, bootstrap_requested=50,
        bootstrap_connected=len(boot), bootstrap_disconnected=disconnected,
        bootstrap_single_block_windows=sum(n == 1 for n in block_counts),
        windows_sampled=len(windows), overlap_components=len(result["mbar_components"]),
        overlap_edge_threshold=result["overlap_edge_threshold"],
        accepted_production_ps_total=sum(block_counts) * 5,
        accepted_production_ps_min=min(block_counts) * 5,
        accepted_production_ps_max=max(block_counts) * 5,
        raw_frames=sum(w["raw_frames"] for w in summaries),
        decorrelated_frames=sum(w["mbar_frames"] for w in summaries),
        decorrelated_frames_min=min(w["mbar_frames"] for w in summaries),
        decorrelated_frames_max=max(w["mbar_frames"] for w in summaries),
        unbiased_target_ESS=float(1 / np.sum(weights**2)),
        changed_owner_frames=sum(w["changed_owner_frames"] for w in summaries),
        block_drift_windows=sum(w["block_drift_flag"] for w in summaries),
        maximum_relative_ES_kcal_mol=float(pmf[maximum]),
        maximum_bin_s_A=float(mid[maximum]),
        maximum_fixed_bin_bootstrap95_kcal_mol=[float(low[maximum]), float(high[maximum])],
        activation_barrier_claim=False, sampling_convergence=False,
        software={p:importlib.metadata.version(p) for p in ["numpy","scipy","pymbar","matplotlib"]},
    )
    (out / "pmf_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    ax.fill_between(mid, low, high, color="#98bbdb", alpha=.65, label="95% block-bootstrap interval")
    ax.plot(mid, pmf, "o-", color="#164b73", markersize=3, label="Native two-CV MBAR marginal")
    ax.set(xlabel="Projected path arc length s (Å); ES → acylenzyme",
           ylabel="F(s) − F(first ES bin) (kcal/mol)", title="Conditional 1D path PMF (25 equal-width bins)")
    ax.axhline(0, color="gray", linewidth=.7)
    ax.legend(fontsize=8)
    fes = result["diagnostic_fes"]
    values = np.array([[np.nan if x is None else x for x in row]
                       for row in fes["free_energy_kcal_mol"]])
    image = axes[1].pcolormesh(fes["x_edges"], fes["y_edges"], values.T, cmap="viridis", shading="flat")
    axes[1].plot(centers[:, 0], centers[:, 1], color="white", linewidth=.8, alpha=.8)
    axes[1].set(xlabel="q_acyl = C−N − C−OG (Å)", ylabel="q_H = OG−H − N−H (Å)", title="Native two-CV free-energy density (minimum = 0)")
    axes[1].set_facecolor("#dddddd")
    fig.colorbar(image, ax=axes[1], label="Relative free energy (kcal/mol)")
    fig.suptitle(f"PRELIMINARY · 300 K DFTB3 · fixed 12 Å exterior · 63/63 windows · target ESS {summary['unbiased_target_ESS']:.1f}", fontsize=11)
    fig.text(.5, .015, "5 ps block resampling; sampling convergence and an activation free energy are not established.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.95))
    for ext in ["png", "svg"]:
        fig.savefig(out / f"conditional_path_pmf.{ext}", dpi=180)
    np.savez_compressed(out / "path_pmf_bootstrap.npz", edges=edges, pmf=pmf, bootstrap=boot, weights=weights, s=coordinate)
    (out / "scientific_acceptance.json").write_text(json.dumps({"process_outputs_complete":True,"connected_sampled_overlap":True,"all63_sampled":len(windows)==63,"activation_barrier_validated":False,"sampling_converged":False,"immutable_inputs":"input_source_manifest.json"},indent=2)+"\n")
    logging.info("Finished requested PMF: %s", summary)


if __name__ == "__main__":
    main()
