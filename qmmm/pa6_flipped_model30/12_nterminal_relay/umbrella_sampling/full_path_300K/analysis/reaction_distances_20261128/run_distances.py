"""Reconstruct direct distances against the exact frozen pooled weight order."""

from __future__ import annotations

import csv
import json
import logging
from collections import defaultdict
from pathlib import Path

import numpy as np

from nyl12_covalent.reaction_distance_stats import DISTANCES, extract_row, sha256
from nyl12_covalent.relay_umbrella_projection import project_path
from plot_distances import plot_all

WORK = Path(__file__).resolve().parent
BASELINE = WORK.parent / "review_20261128_01/checkpoint_001_main70ps_supplement43"
LOGGER = logging.getLogger(__name__)
OWNERS = ("hg1_owner", "nalpha_h1_owner", "nalpha_h2_owner")
EXTRA = ("r_nalpha_h1", "r_nalpha_h2")


def save_json(name: str, value: object) -> None:
    """Save finite derived evidence."""
    (WORK / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def summarize(values: np.ndarray, weights: np.ndarray) -> dict:
    """Return empirical quantiles and weighted population spread, never a CI."""
    order = np.argsort(values, kind="stable")
    cdf = np.cumsum(weights[order]) / weights.sum()
    quantiles = values[order][np.searchsorted(cdf, [0.1, 0.5, 0.9])]
    return dict(
        mean_A=float(np.average(values, weights=weights)),
        p10_A=float(quantiles[0]),
        median_A=float(quantiles[1]),
        p90_A=float(quantiles[2]),
    )


def main() -> None:
    """Verify all input hashes and emit ordered physical-coordinate diagnostics."""
    sources = json.loads((BASELINE / "sources.json").read_text())
    weights = np.load(BASELINE / "weights.npz")
    proofs = defaultdict(list)
    for proof in sources["proofs"]:
        proofs[proof["state_id"]].append(proof)
    hashes = {
        str(BASELINE / n): sha256(BASELINE / n) for n in ("sources.json", "weights.npz")
    }
    rows, native, state_offsets = [], [], []
    for state_index, state in enumerate(sources["states"]):
        block_rows = []
        for proof in proofs[state["state_id"]]:
            path = Path(
                next(
                    p for p in proof["source_sha256"] if p.endswith("/raw_features.csv")
                )
            )
            digest = sha256(path)
            if digest != proof["source_sha256"][str(path)]:
                raise ValueError(f"Source hash mismatch: {path}")
            hashes[str(path)] = digest
            with path.open() as handle:
                raw_rows = list(csv.DictReader(handle))
            if len(raw_rows) != 100:
                raise ValueError(f"Expected 100 raw rows: {path}")
            for raw in raw_rows:
                row = extract_row(raw)
                row.update({k: float(raw[k]) for k in EXTRA})
                row.update(
                    {k: float(raw[k]) for k in ("q_acyl_A", "q_act_A", "q_del_A")}
                )
                row.update(
                    state_index=state_index,
                    state_id=state["state_id"],
                    source=str(path),
                    source_frame=int(raw["frame"]),
                )
                block_rows.append(row)
        times = np.array([r["time_ps"] for r in block_rows])
        if not np.allclose(np.diff(times), 0.05, atol=1e-4, rtol=0):
            raise ValueError(f"Noncontiguous blocks: {state['state_id']}")
        origin = times[-1] - state["production_horizon_ps"]
        for row in block_rows:
            row["production_time_ps"] = round(row["time_ps"] - origin, 4)
        selected = state["window_statistics"]["selected_indices"]
        if selected != np.flatnonzero(times - origin > 5 + 1e-5).tolist():
            raise ValueError("Selection is not exactly one five-ps discard")
        use = [block_rows[i] for i in selected]
        if len(use) != weights["counts"][state_index]:
            raise ValueError("Frozen state count mismatch")
        state_offsets.append(len(rows))
        rows.extend(use)
        native.extend(block_rows)
        LOGGER.info("Verified %s: %d retained rows", state["state_id"], len(use))
    q = np.array([[r[k] for k in ("q_acyl_A", "q_act_A", "q_del_A")] for r in rows])
    np.testing.assert_allclose(q, weights["q"], atol=1e-12, rtol=0)
    projection = project_path(q, np.asarray(sources["config"]["path_3cv_A"]))
    s = projection["s"]
    edges = np.array(sources["config"]["s_bin_edges"])
    bins = np.searchsorted(edges, s, side="right") - 1
    centers = (edges[:-1] + edges[1:]) / 2
    lw = weights["logweights"]
    data = {k: np.array([r[k] for r in rows]) for k in (*DISTANCES, *EXTRA, *OWNERS)}
    diagnostics, distances = [], []
    for b in range(len(centers)):
        mask = bins == b
        if not mask.any():
            continue
        w = np.exp(lw[mask] - lw[mask].max())
        local_rows = np.flatnonzero(mask)
        dominant = int(local_rows[np.argmax(w)])
        item = dict(
            s=float(centers[b]),
            n=int(mask.sum()),
            weight_ESS=float(w.sum() ** 2 / (w @ w)),
            dominant_frame_index=dominant,
            dominant_fraction=float(w.max() / w.sum()),
            dominant_row=rows[dominant],
            transverse_A=summarize(projection["transverse_A"][mask], w),
            owners={},
        )
        for key in OWNERS:
            item["owners"][key] = {
                str(int(owner)): dict(
                    weighted=float(w[data[key][mask] == owner].sum() / w.sum()),
                    unweighted=float(np.mean(data[key][mask] == owner)),
                )
                for owner in np.unique(data[key][mask])
            }
        diagnostics.append(item)
        for key in (*DISTANCES, *EXTRA):
            distances.append(
                dict(
                    s=float(centers[b]),
                    distance=key,
                    weighted=summarize(data[key][mask], w),
                    unweighted=summarize(data[key][mask], np.ones(len(w))),
                )
            )
    save_json("bin_diagnostics.json", diagnostics)
    save_json("distance_statistics.json", distances)
    save_json("source_hashes.json", hashes)
    with (WORK / "retained_rows.csv").open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]) + ["s", "transverse_A", "logweight"]
        )
        writer.writeheader()
        for i, row in enumerate(rows):
            writer.writerow(
                dict(
                    row,
                    s=s[i],
                    transverse_A=projection["transverse_A"][i],
                    logweight=lw[i],
                )
            )
    plot_all(WORK, rows, native, data, s, centers, distances)
    save_json(
        "scientific_acceptance.json",
        dict(
            frozen_baseline=str(BASELINE),
            rows=len(rows),
            states=len(sources["states"]),
            direct_distances=5,
            CV_identity_checks=3,
            exact_q_order_verified=True,
            source_hashes_verified=len(hashes),
            discard_once_ps=5,
            intervals="Weighted empirical 10–90% distributions; not confidence intervals",
            limitations="Biased histories remain correlated; sampling convergence unresolved. "
            "Distance geometry alone does not establish bond order, product identity or a TS.",
        ),
    )
    LOGGER.info("Completed exact pooled diagnostics and plots")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
