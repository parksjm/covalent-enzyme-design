"""Uncensored completed-block collection and bias-only two-dimensional MBAR."""

from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import logsumexp

from .neb_umbrella_bias import HarmonicWindow, bias_energy

KB = 0.00198720425864083


def digest(path: Path) -> str:
    """Hash a compact source or observable artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collect(prepared: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Read only explicitly declared, completed production blocks; censor no chemistry."""
    cfg = json.loads(prepared.read_text())
    if digest(prepared) != (prepared.parent / "prepared.sha256").read_text().strip():
        raise ValueError("Prepared campaign identity changed")
    windows = []
    for window in cfg["windows"]:
        blocks = []
        hashes = {}
        for stage in window["stages"]:
            if stage["kind"] != "production":
                continue
            folder = prepared.parent / window["label"] / stage["label"]
            paths = [
                folder / "acceptance.json",
                folder / "samples_and_cross_bias.npz",
                folder / "observables/observables.json",
            ]
            if not all(p.is_file() for p in paths):
                continue
            accepted = json.loads(paths[0].read_text())
            if (
                not accepted["finite_energy_temperature"]
                or accepted["completed_steps"] <= 0
            ):
                raise ValueError(f"Invalid completed production acceptance: {folder}")
            receipt = json.loads(paths[2].read_text())
            csv_path = folder / "observables/observables.csv"
            if (
                digest(csv_path) != receipt["csv_sha256"]
                or receipt["frames_censored"] != 0
            ):
                raise ValueError("Observable identity/censoring contract changed")
            with csv_path.open() as stream:
                rows = list(csv.DictReader(stream))
            if len(rows) != receipt["frames"] or not rows:
                raise ValueError("Observable frame count mismatch")
            samples = np.array(
                [
                    [
                        float(r[k])
                        for k in ("q_acyl", "q_h", "cn_progress", "cn_nalpha_hpool")
                    ]
                    for r in rows
                ]
            )
            if not np.isfinite(samples).all():
                raise ValueError("Nonfinite production observables")
            with np.load(paths[1]) as archive:
                if not np.allclose(
                    archive["cv_values"], samples[:, :2], rtol=0, atol=1e-8
                ):
                    raise ValueError("Cross-bias/observable frame identity mismatch")
            blocks.append(dict(label=stage["label"], samples=samples, rows=rows))
            hashes.update({str(p): digest(p) for p in [*paths, csv_path]})
        if blocks:
            windows.append(
                dict(
                    label=window["label"],
                    bias=window["bias"],
                    blocks=blocks,
                    source_sha256=hashes,
                )
            )
    return cfg, windows


def components(matrix: np.ndarray, threshold: float = 1e-3) -> list[list[int]]:
    """Return undirected overlap components at an explicit diagnostic edge threshold."""
    graph = (matrix > threshold) | (matrix.T > threshold)
    remaining = set(range(len(matrix)))
    result = []
    while remaining:
        group = {remaining.pop()}
        front = list(group)
        while front:
            linked = set(np.flatnonzero(graph[front.pop()])) & remaining
            remaining -= linked
            group |= linked
            front.extend(linked)
        result.append(sorted(map(int, group)))
    return result


def histogram_fes(
    points: np.ndarray, weights: np.ndarray, temperature: float, bins: int
) -> dict[str, Any]:
    """Return weighted histogram density/FES with unsampled bins explicitly missing."""
    hist, x, y = np.histogram2d(points[:, 0], points[:, 1], bins=bins, weights=weights)
    density = hist / np.outer(np.diff(x), np.diff(y))
    populated = density > 0
    fes = np.full(hist.shape, np.nan)
    fes[populated] = -KB * temperature * np.log(density[populated])
    fes[populated] -= fes[populated].min()
    return dict(
        x_edges=x.tolist(),
        y_edges=y.tolist(),
        probability=hist.tolist(),
        free_energy_kcal_mol=[
            [float(v) if np.isfinite(v) else None for v in row] for row in fes
        ],
        unsampled_bins=int((~populated).sum()),
        normalization="Histogram probability divided by bin area; minimum populated density-derived F set0",
        barrier_claim=False,
    )


def analyze(
    cfg: dict[str, Any], windows: list[dict[str, Any]], bins: int = 20
) -> dict[str, Any]:
    """Estimate interim unbiased weights under the explicitly shared physical Hamiltonian.

    All valid rows enter raw coverage/state diagnostics. MBAR uses deterministic
    per-window statistical thinning, never chemical-state selection. U0 cancels
    between states, so sampled reduced potentials are beta*bias; target u=0.
    """
    from pymbar import MBAR, timeseries

    if bins < 2 or not np.isfinite(cfg["temperature_K"]) or cfg["temperature_K"] <= 0:
        raise ValueError("Positive temperature and at least2bins required")
    result = dict(
        status="no_completed_production" if not windows else "interim_incomplete",
        windows_total=len(cfg["windows"]),
        windows_sampled=len(windows),
        sampling_convergence=False,
        activation_barrier_claim=False,
        chemical_frames_censored=0,
        window_summaries=[],
        hysteresis="Not established: independent forward/reverse replica sampling not provided",
    )
    if not windows:
        return result
    raw = [np.concatenate([b["samples"] for b in w["blocks"]]) for w in windows]
    selected = []
    for window, samples in zip(windows, raw, strict=True):
        g = [
            1.0
            if len(samples) < 4 or np.ptp(samples[:, i]) == 0
            else float(timeseries.statistical_inefficiency(samples[:, i]))
            for i in range(3)
        ]
        if not np.isfinite(g).all() or min(g) < 1:
            raise ValueError("Invalid statistical inefficiency")
        stride = int(np.ceil(max(g)))
        chosen = samples[::stride]
        selected.append(chosen)
        owners = Counter(
            (r["hg1_owner"], r["nalpha_h1_owner"], r["nalpha_h2_owner"])
            for b in window["blocks"]
            for r in b["rows"]
        )
        means = [b["samples"][:, :3].mean(axis=0) for b in window["blocks"]]
        drift = bool(
            len(means) > 1
            and np.any(np.abs(means[-1] - means[0]) > samples[:, :3].std(axis=0))
        )
        result["window_summaries"].append(
            dict(
                label=window["label"],
                blocks=len(window["blocks"]),
                raw_frames=len(samples),
                g_q_acyl=g[0],
                g_q_h=g[1],
                g_proton_PI=g[2],
                decorrelation_stride=stride,
                mbar_frames=len(chosen),
                thinning_policy="Fixed stride from maxg over qAcyl,qH,protonPI; all chemical states eligible",
                owner_states=[
                    dict(owners=list(k), frames=v) for k, v in owners.items()
                ],
                changed_owner_frames=sum(
                    r["owner_changed"].lower() == "true"
                    for b in window["blocks"]
                    for r in b["rows"]
                ),
                block_drift_flag=drift,
                proton_PI_negative_frames=int((samples[:, 2] < 0).sum()),
                proton_PI_positive_frames=int((samples[:, 2] > 0).sum()),
                nalpha_CN_range=[
                    float(samples[:, 3].min()),
                    float(samples[:, 3].max()),
                ],
            )
        )
    points = np.concatenate(raw)
    edges = [
        np.linspace(points[:, i].min() - 0.001, points[:, i].max() + 0.001, bins + 1)
        for i in range(2)
    ]
    hist = [
        np.histogram2d(s[:, 0], s[:, 1], bins=edges)[0].ravel() / len(s) for s in raw
    ]
    root_hist = np.sqrt(hist)
    empirical = root_hist @ root_hist.T
    result.update(
        coverage_points=points[:, :3].tolist(),
        histogram_overlap=empirical.tolist(),
        histogram_components=components(empirical),
        source_sha256={k: v for w in windows for k, v in w["source_sha256"].items()},
    )
    if len(windows) < 2 or any(len(s) < 2 for s in selected):
        result["mbar_status"] = "insufficient_sampled_windows_or_decorrelated_frames"
        return result
    sampled_counts = {
        w["label"]: len(s) for w, s in zip(windows, selected, strict=True)
    }
    counts = np.array([sampled_counts.get(w["label"], 0) for w in cfg["windows"]])
    samples = np.concatenate(selected)
    beta = 1 / (KB * cfg["temperature_K"])
    all_biases = [HarmonicWindow(**w["bias"]) for w in cfg["windows"]]
    reduced = np.array(
        [[beta * bias_energy(tuple(p[:2]), w) for p in samples] for w in all_biases]
    )
    mbar = MBAR(reduced, counts, verbose=False, solver_protocol="robust")
    logw = -logsumexp(
        np.log(counts[counts > 0])[:, None]
        + mbar.f_k[counts > 0, None]
        - reduced[counts > 0],
        axis=0,
    )
    weights = np.exp(logw - logsumexp(logw))
    overlap = mbar.compute_overlap()["matrix"]
    result.update(
        mbar_status="interim_bias_only_common_U0",
        mbar_window_labels=[w["label"] for w in cfg["windows"]],
        mbar_N_k=counts.tolist(),
        mbar_overlap=overlap.tolist(),
        mbar_components=components(overlap),
        connected=len(components(overlap[np.ix_(counts > 0, counts > 0)])) == 1,
        connected_scope="Sampled windows only; unsampled windows remain explicit N_k=0",
        overlap_edge_threshold=1e-3,
        unbiased_effective_sample_size=float(1 / np.sum(weights**2)),
        target_reduced_potential=0,
        common_U0_cancellation="Same topology/QM mask/charge/temperature/frozen boundary; only exact piecewise Amber biases differ",
        primary_axes=["q_acyl_A", "proton_PI_coordination_difference"],
        diagnostic_axes=["q_acyl_A", "q_h_A"],
    )
    if not result["connected"]:
        result.update(
            status="interim_disconnected_overlap",
            fes_status="suppressed_disconnected_sampled_overlap",
            relative_basin_claim=False,
            interpretation="Disconnected sampled-state overlap: only coverage/overlap/state diagnostics; no global FES or relative basin comparison.",
        )
        return result
    result.update(
        fes_status="interim_connected_overlap_not_convergence",
        primary_fes=histogram_fes(
            samples[:, [0, 2]], weights, cfg["temperature_K"], bins
        ),
        diagnostic_fes=histogram_fes(
            samples[:, :2], weights, cfg["temperature_K"], bins
        ),
    )
    if len(windows) == len(cfg["windows"]):
        result["status"] = "interim_all_windows_sampled_not_converged"
    return result
