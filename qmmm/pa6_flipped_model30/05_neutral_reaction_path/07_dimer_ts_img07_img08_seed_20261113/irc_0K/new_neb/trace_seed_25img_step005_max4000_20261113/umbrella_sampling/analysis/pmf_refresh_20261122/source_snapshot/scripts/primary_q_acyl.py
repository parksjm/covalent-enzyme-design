"""Primary conditional q_acyl marginal from unchanged native two-CV MBAR."""
from __future__ import annotations
import csv
import hashlib
import json
import logging
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from requested_path_pmf import KB, HarmonicWindow, collect, solve


def density_profile(x: np.ndarray, weights: np.ndarray, edges: np.ndarray,
                    beta: float, reference: int) -> np.ndarray:
    """Accumulate positive bin masses directly without cumulative subtraction."""
    if np.any(x < edges[0]) or np.any(x > edges[-1]):
        raise ValueError("Coordinate outside declared bounds")
    indices = np.minimum(np.searchsorted(edges, x, side="right") - 1, len(edges) - 2)
    mass = np.bincount(indices, weights=weights, minlength=len(edges) - 1)
    if mass[reference] <= 0:
        raise ValueError("Empty fixed ES reference bin")
    result = np.full(len(mass), np.nan)
    valid = mass > 0
    result[valid] = -np.log(mass[valid] / np.diff(edges)[valid]) / beta
    return result - result[reference]


def regression() -> None:
    """Verify fifty-decade probability range, sample permutation and ES zero."""
    x = np.array([.5, 1.5, 2.5])
    w = np.array([1., 1e-30, 1e-50])
    edges = np.arange(4.)
    expected = -np.log(w)
    result = density_profile(x, w, edges, 1., 0)
    np.testing.assert_allclose(result, expected)
    np.testing.assert_allclose(density_profile(x[::-1], w[::-1], edges, 1., 0), expected)
    assert np.isfinite(result).all()


def main() -> None:
    """Refit completed-block bootstrap and write primary natural-coordinate PMF."""
    logging.basicConfig(level=logging.INFO)
    regression()
    p = Path(__file__).resolve().parents[2]
    cfg, windows = collect(p / "frozen_inputs/prepared.json")
    result = json.loads((p / "mbar/latest.json").read_text())
    summaries = result["window_summaries"]
    strides = [w["decorrelation_stride"] for w in summaries]
    raw = [np.concatenate([b["samples"] for b in w["blocks"]]) for w in windows]
    selected = [a[::s] for a, s in zip(raw, strides, strict=True)]
    beta = 1 / (KB * cfg["temperature_K"])
    biases = [HarmonicWindow(**w["bias"]) for w in cfg["windows"]]
    weights, free, connected = solve(selected, biases, beta)
    assert connected
    all_x = np.concatenate(raw)[:, 0]
    previous = np.load(p / "frozen_inputs/previous_profile.npz")
    original_edges = previous["edges"]
    width = original_edges[1] - original_edges[0]
    first = min(0, int(np.floor((all_x.min() - original_edges[0]) / width)))
    last = max(len(original_edges) - 1, int(np.ceil((all_x.max() - original_edges[0]) / width)))
    edges = original_edges[0] + np.arange(first, last + 1) * width
    reference = int(np.searchsorted(edges, cfg["windows"][0]["q_acyl_center"], side="right") - 1)
    x = np.concatenate(selected)[:, 0]
    pmf = density_profile(x, weights, edges, beta, reference)
    rng = np.random.default_rng(20261121)
    bootstrap = []
    disconnected = 0
    for replicate in range(50):
        resampled = []
        for window, stride in zip(windows, strides, strict=True):
            blocks = window["blocks"]
            picks = rng.integers(0, len(blocks), size=len(blocks))
            resampled.append(np.concatenate([blocks[i]["samples"] for i in picks])[::stride])
        bw, _, ok = solve(resampled, biases, beta, free)
        if not ok:
            disconnected += 1
            continue
        bootstrap.append(density_profile(np.concatenate(resampled)[:, 0], bw, edges, beta, reference))
        logging.info("Primary q_acyl bootstrap %d/50", replicate + 1)
    boot = np.array(bootstrap)
    low, high = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    mid = (edges[:-1] + edges[1:]) / 2
    maximum = int(np.nanargmax(pmf))
    bcounts = [len(w["blocks"]) for w in windows]
    summary = dict(
        coordinate="q_acyl=d(C,N_leaving)-d(C,OG1), Angstrom; native two-CV MBAR marginal over q_H and all remaining sampled coordinates",
        reference="Density of fixed q_acyl histogram bin containing first ES umbrella center; not an integrated ES basin",
        reference_bin_edges_A=[float(edges[reference]),float(edges[reference+1])],
        bins=len(edges)-1, temperature_K=cfg["temperature_K"], Hamiltonian="DFTB3/3ob, frozen original12A exterior, historical neutral-Thr-N ASP272/ASP274 declared proton inventory",
        uncertainty="50 within-window5ps block-bootstrap refits; fixed original statistical thinning; pointwise95percentile interval; no convergence claim",
        bootstrap_connected=len(boot),bootstrap_disconnected=disconnected,
        windows_sampled=len(windows),windows_total=len(cfg["windows"]),
        overlap_components=len(result["mbar_components"]),overlap_edge_threshold=.001,
        accepted_production_ps_total=5*sum(bcounts),accepted_production_ps_min=5*min(bcounts),accepted_production_ps_max=5*max(bcounts),
        raw_frames=sum(len(s) for s in raw),decorrelated_frames=len(weights),
        unbiased_target_ESS=float(1/np.sum(weights**2)),
        block_drift_windows=[s["label"] for s in summaries if s["block_drift_flag"]],
        maximum_relative_ES_bin_kcal_mol=float(pmf[maximum]),
        maximum_q_acyl_A=float(mid[maximum]),
        fixed_maximum_bin_bootstrap95_kcal_mol=[float(low[maximum]),float(high[maximum])],
        unsampled_bins=int(np.isnan(pmf).sum()),
        sampling_convergence=False,activation_barrier_claim=False,
        numerical_accumulation="Direct positive np.bincount masses; regression retains weights1,1e-30,1e-50 and sample-order invariance",
        input_source_manifest_sha256=hashlib.sha256((p/"input_source_manifest.json").read_bytes()).hexdigest())
    (p / "primary_q_acyl_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    with (p/"primary_q_acyl_pmf.csv").open("w") as stream:
        writer=csv.writer(stream);writer.writerow(["q_acyl_A","F_relative_ES_bin_kcal_mol","bootstrap95_low","bootstrap95_high"]);writer.writerows(zip(mid,pmf,low,high,strict=True))
    np.savez_compressed(p/"primary_q_acyl_bootstrap.npz",edges=edges,pmf=pmf,bootstrap=boot)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    ax=axes[0]
    ax.fill_between(mid,low,high,color="#98bbdb",alpha=.65,label="95% block-bootstrap interval")
    ax.plot(mid,pmf,"o-",color="#164b73",markersize=3,label="Native two-CV MBAR marginal")
    ax.axhline(0,color="gray",linewidth=.7)
    ax.axvspan(edges[reference],edges[reference+1],color="#d4e8d1",alpha=.45,label="Fixed ES reference bin")
    ax.set(xlabel="q_acyl = d(C,N) − d(C,OG) (Å)",ylabel="F(q_acyl) − F(ES reference bin) (kcal/mol)",title="Conditional 1D PMF (unchanged bin width and ES reference)")
    ax.legend(fontsize=8,loc="upper left")
    ax.text(.98, .96, f"Apparent maximum: {pmf[maximum]:.1f} kcal/mol", transform=ax.transAxes, ha="right", va="top", fontsize=9)
    fes=result["diagnostic_fes"]
    values=np.array([[np.nan if v is None else v for v in row] for row in fes["free_energy_kcal_mol"]])
    image=axes[1].pcolormesh(fes["x_edges"],fes["y_edges"],values.T,cmap="viridis",shading="flat")
    axes[1].set(xlabel="q_acyl = C−N − C−OG (Å)",ylabel="q_H = OG−H − N−H (Å)",title="Native two-CV free-energy density (minimum = 0)",facecolor="#dddddd")
    fig.colorbar(image,ax=axes[1],label="Relative free energy (kcal/mol)")
    fig.suptitle(f"PRELIMINARY · 300 K DFTB3 · fixed 12 Å exterior · 63/63 windows · target ESS {summary['unbiased_target_ESS']:.1f}",fontsize=11)
    fig.text(.5,.015,"95% intervals from 50 refitted 5 ps block resamples; activation free energy and convergence remain unestablished.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.95))
    for ext in ["png","svg"]:fig.savefig(p/f"primary_q_acyl_pmf.{ext}",dpi=180)
    (p/"primary_scientific_acceptance.json").write_text(json.dumps({"process_outputs_complete":True,"numerical_dynamic_range_regression_passed":True,"connected_sampled_overlap":connected,"all63_sampled":len(windows)==63,"activation_barrier_validated":False,"sampling_converged":False},indent=2)+"\n")
    logging.info("Primary result %s",summary)


if __name__ == "__main__":
    main()
