"""Diagnose unchanged pooled weights: corners, support, lineage and coordinate."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import logging
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import logsumexp

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K"
)
SOURCE = BASE / "analysis/review_20261128_01/current"
OUT = Path(__file__).parent
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.relay_umbrella_projection import project_path  # noqa: E402

HELPER = (
    BASE
    / "analysis/fixed_burnin_20261126/projection_sensitivity_20261126/projection.py"
)
SPEC = importlib.util.spec_from_file_location("frozen_projection", HELPER)
assert SPEC and SPEC.loader
PROJECTION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PROJECTION)
KT = 0.00198720425864083 * 300


def save(name: str, rows: list[dict[str, Any]]) -> None:
    """Write diagnostic rows as JSON authority and CSV companion."""
    (OUT / f"{name}.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False) + "\n"
    )
    if rows:
        with (OUT / f"{name}.csv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main() -> None:
    """Evaluate exact displayed samples without refitting their weights."""
    metadata = json.loads((SOURCE / "sources.json").read_text())
    point = json.loads((SOURCE / "point.json").read_text())
    data = np.load(SOURCE / "weights.npz")
    q, lw, counts = data["q"], data["logweights"], data["counts"]
    states = metadata["states"]
    owner = np.repeat(np.arange(len(counts)), counts)
    age = np.concatenate(
        [
            np.asarray(s["window_statistics"]["selected_indices"])
            * s["window_statistics"]["frame_interval_ps"]
            + s["window_statistics"]["selected_start_ps"]
            - s["window_statistics"]["selected_indices"][0]
            * s["window_statistics"]["frame_interval_ps"]
            for s in states
        ]
    )
    assert len(q) == 82000 and len(states) == 100 and len(age) == len(q)
    path = np.asarray(metadata["config"]["path_3cv_A"])
    hard = project_path(q, path)
    s = hard["s"]
    edges = np.asarray(metadata["config"]["s_bin_edges"])
    centers = (edges[1:] + edges[:-1]) / 2
    bins = np.searchsorted(edges, s, side="right") - 1
    corners = np.r_[0, np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
    length = corners[-1]
    corners /= length
    corner_id = np.argmin(abs(s[:, None] - corners), axis=1)
    on_corner = abs(s - corners[corner_id]) < 1e-12
    # Exclude the endpoints, whose terminal segments extend beyond the path.
    on_corner &= (corner_id > 0) & (corner_id < len(path) - 1)
    w = np.exp(lw - logsumexp(lw))
    rows, contributions, dominant = [], [], []
    corner_rows = []
    for j, corner in enumerate(corners):
        if not -0.05 <= corner <= 0.3:
            continue
        use = on_corner & (corner_id == j)
        surrounding = bins == np.searchsorted(edges, corner, side="right") - 1
        corner_rows.append(
            dict(
                vertex=j,
                s=float(corner),
                raw_count=int(use.sum()),
                probability_mass=float(w[use].sum()),
                fraction_of_bin_mass=float(w[use].sum() / w[surrounding].sum()),
                q_acyl=float(path[j, 0]),
                q_act=float(path[j, 1]),
                q_del=float(path[j, 2]),
            )
        )
    for i, center in enumerate(centers):
        if not -0.05 <= center <= 0.3:
            continue
        use = bins == i
        ids = np.flatnonzero(use)
        local = w[use] / w[use].sum()
        group_mass: dict[str, float] = {}
        for j in np.unique(owner[use]):
            mass = float(w[use & (owner == j)].sum() / w[use].sum())
            st = states[j]
            contributions.append(
                dict(
                    s=float(center),
                    state=st["state_id"],
                    group=st["group"],
                    mass_fraction=mass,
                )
            )
            group_mass[st["group"]] = group_mass.get(st["group"], 0) + mass
        ranked = ids[np.argsort(w[ids])[-5:][::-1]]
        for rank, idx in enumerate(ranked, 1):
            dominant.append(
                dict(
                    s=float(center),
                    rank=rank,
                    frame_index=int(idx),
                    state=states[owner[idx]]["state_id"],
                    age_ps=float(age[idx]),
                    weight_fraction=float(w[idx] / w[use].sum()),
                    q_acyl=float(q[idx, 0]),
                    q_act=float(q[idx, 1]),
                    q_del=float(q[idx, 2]),
                    projected_s=float(s[idx]),
                    segment=int(hard["segment"][idx]),
                    on_vertex=bool(on_corner[idx]),
                )
            )
        rows.append(
            dict(
                s=float(center),
                left=float(edges[i]),
                right=float(edges[i + 1]),
                pmf_kcal_mol=point["projections"]["hard"]["pmf"][i],
                raw_count=int(use.sum()),
                kish_ess=float(1 / (local @ local)),
                max_frame_fraction=float(local.max()),
                vertex_raw_count=int((use & on_corner).sum()),
                vertex_mass_fraction=float(w[use & on_corner].sum() / w[use].sum()),
                main57_mass_fraction=group_mass.get("main57", 0),
                mean_age_ps=float(local @ age[use]),
                age_le10ps_fraction=float(w[use & (age <= 10)].sum() / w[use].sum()),
                age_gt35ps_fraction=float(w[use & (age > 35)].sum() / w[use].sum()),
                dominant_state=states[owner[ranked[0]]]["state_id"],
            )
        )
    save("early_bins", rows)
    save("early_vertices", corner_rows)
    save("lineage_contributions", contributions)
    save("dominant_frames", dominant)
    nodes, _ = PROJECTION.uniform_nodes(path, 153)
    bgp, gradient = PROJECTION.smooth_progress(q, nodes, 56.36835607868901)
    gnorm = np.linalg.norm(gradient, axis=1)
    bgp_edges = np.linspace(0, 1, 41)
    hb = PROJECTION.local_histogram(bgp, lw, bgp_edges)
    actual = -KT * (np.asarray(hb["logdensity"], float) - hb["logdensity"][0])
    np.testing.assert_allclose(
        actual, np.asarray(point["projections"]["BGP"]["pmf"], float), equal_nan=True
    )
    np.testing.assert_array_equal(
        np.bincount(bins, minlength=len(centers)),
        point["projections"]["hard"]["counts"],
    )
    bgp_rows = []
    for i, center in enumerate((bgp_edges[1:] + bgp_edges[:-1]) / 2):
        if center > 0.3:
            continue
        use = (bgp >= bgp_edges[i]) & (bgp < bgp_edges[i + 1])
        local = w[use] / w[use].sum()
        bgp_rows.append(
            dict(
                s=float(center),
                pmf_kcal_mol=float(actual[i]),
                kish_ess=float(hb["ess"][i]),
                raw_count=int(use.sum()),
                gradient_norm_weighted_mean=float(local @ gnorm[use]),
                gradient_norm_min=float(gnorm[use].min()),
                gradient_norm_max=float(gnorm[use].max()),
            )
        )
    save("bgp_gradient", bgp_rows)
    qa_edges = np.arange(-2.5, 3.00001, 0.1)
    qa_centers = (qa_edges[1:] + qa_edges[:-1]) / 2
    physical = point["projections"]["initial_q_act"]
    physical_rows = [
        dict(
            q_act=float(x),
            pmf_kcal_mol=physical["pmf"][i],
            raw_count=physical["counts"][i],
            kish_ess=physical["ESS"][i],
            max_frame_fraction=physical["maximum_frame_weight_fraction"][i],
        )
        for i, x in enumerate(qa_centers)
        if -1.5 <= x <= 0.8
    ]
    save("initial_branch_qact", physical_rows)
    plt.rcParams.update(
        {"font.size": 16, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, ax = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    x = [r["s"] for r in rows]
    ax[0, 0].plot(
        x, [r["pmf_kcal_mol"] for r in rows], "o-", label="Hard closest segment"
    )
    ax[0, 0].plot(
        [r["s"] for r in bgp_rows],
        [r["pmf_kcal_mol"] for r in bgp_rows],
        "o-",
        label="BGP; different density measure",
    )
    for r in corner_rows:
        if r["raw_count"]:
            ax[0, 0].axvline(r["s"], color=".65", alpha=0.5, lw=0.8)
    ax[0, 0].set(xlabel="Named path coordinate", ylabel="Own-reference PMF (kcal/mol)")
    ax[0, 0].legend(fontsize=13)
    ax[0, 1].plot(
        x, [r["vertex_mass_fraction"] for r in rows], "o-", label="At exact vertices"
    )
    ax[0, 1].plot(
        x, [r["max_frame_fraction"] for r in rows], "s-", label="Largest frame"
    )
    ax[0, 1].set(
        xlabel="Hard path coordinate s",
        ylabel="Fraction of bin probability",
        ylim=(0, 1),
    )
    ax[0, 1].legend(fontsize=13)
    ax[1, 0].semilogy(x, [r["raw_count"] for r in rows], "o-", label="Raw frames")
    ax[1, 0].semilogy(x, [r["kish_ess"] for r in rows], "s-", label="Kish weight ESS")
    ax[1, 0].set(xlabel="Hard path coordinate s", ylabel="Count / weight ESS")
    ax[1, 0].legend(fontsize=13)
    ax[1, 1].plot(qa_centers, np.asarray(physical["pmf"], float), "o-")
    ax[1, 1].set(
        xlim=(-1.5, 0.8),
        ylim=(-1, 17),
        xlabel="Physical q_act (Å)",
        ylabel="Conditional PMF (kcal/mol)",
        title="Initial branch: −2 ≤ q_acyl < −0.8 Å",
    )
    fig.suptitle(
        "Same 82,000 samples and fixed MBAR weights; no smoothing/refit\nCorner atoms, sparse reweighting and coordinate dependence are distinct"
    )
    fig.savefig(OUT / "early_diagnosis.png", dpi=180)
    fig.savefig(OUT / "early_diagnosis.pdf")
    plt.close(fig)
    inputs = [
        SOURCE / "sources.json",
        SOURCE / "weights.npz",
        SOURCE / "point.json",
        HELPER,
        Path(__file__),
        ROOT / "src/nyl12_covalent/relay_umbrella_projection.py",
    ]
    manifest = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    (OUT / "input_sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (OUT / "acceptance.json").write_text(
        json.dumps(
            dict(
                process="complete",
                baseline_counts_match=True,
                baseline_BGP_match=True,
                weights_refit=False,
                sampling_converged=False,
                age_definition="Selected sample age from history selected_indices/start/frame interval; not independent-start age",
                scientific_scope="Coordinate/support diagnostic only; Kish ESS excludes temporal correlation; no barrier acceptance",
            ),
            indent=2,
        )
        + "\n"
    )
    logging.info(
        "Completed early diagnosis: %d bins, %d vertices", len(rows), len(corner_rows)
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
