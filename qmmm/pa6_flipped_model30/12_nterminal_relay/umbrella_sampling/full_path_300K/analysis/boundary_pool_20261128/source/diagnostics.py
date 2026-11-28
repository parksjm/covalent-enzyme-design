"""Uncensored geometry, owner and transverse summaries in exact weight order."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from nyl12_covalent.reaction_distance_stats import DISTANCES, extract_row
from nyl12_covalent.relay_umbrella_projection import project_path
from nyl12_covalent.relay_umbrella_analysis_io import digest


def summarize(data: dict, logweights: np.ndarray) -> dict:
    """Check all three direct-distance CV identities and report both initializations.

    Args:
        data: Source-bound exact106-history pool.
        logweights: Matched unbiased logweights from the unchanged global MBAR fit.

    Returns:
        Per-history descriptive distributions and weighted support diagnostics.
    """
    proofs = defaultdict(list)
    for proof in data["proofs"]:
        proofs[proof["state_id"]].append(proof)
    rows, groups = [], []
    for state in data["states"]:
        native = []
        for proof in proofs[state["state_id"]]:
            path = Path(
                next(
                    k for k in proof["source_sha256"] if k.endswith("/raw_features.csv")
                )
            )
            if digest(path) != proof["source_sha256"][str(path)]:
                raise ValueError("Raw feature identity changed")
            with path.open() as stream:
                for row in csv.DictReader(stream):
                    direct = extract_row(row)
                    direct["q"] = [
                        float(row[k]) for k in ("q_acyl_A", "q_act_A", "q_del_A")
                    ]
                    native.append(direct)
        use = [native[i] for i in state["window_statistics"]["selected_indices"]]
        if len(use) != state["retained_frames"]:
            raise ValueError("Diagnostic row count mismatch")
        rows.extend(use)
        groups.append(
            dict(
                state_id=state["state_id"],
                bias=state["bias"],
                start=state["source_lineage"],
                distances={
                    k: dict(
                        mean_A=float(np.mean([r[k] for r in use])),
                        p10_p50_p90_A=np.percentile(
                            [r[k] for r in use], [10, 50, 90]
                        ).tolist(),
                    )
                    for k in DISTANCES
                },
                owners={
                    k: dict(Counter(str(r[k]) for r in use))
                    for k in ("hg1_owner", "nalpha_h1_owner", "nalpha_h2_owner")
                },
            )
        )
    q = np.array([r["q"] for r in rows])
    np.testing.assert_allclose(q, data["q"], atol=1e-12, rtol=0)
    projection = project_path(q, np.array(data["config"]["path_3cv_A"]))
    edges = np.array(data["config"]["s_bin_edges"])
    bins = np.searchsorted(edges, projection["s"], side="right") - 1
    support = []
    for index in range(len(edges) - 1):
        mask = bins == index
        if not mask.any():
            continue
        w = np.exp(logweights[mask] - logweights[mask].max())
        support.append(
            dict(
                s=float((edges[index] + edges[index + 1]) / 2),
                n=int(mask.sum()),
                weight_Kish_ESS=float(w.sum() ** 2 / (w @ w)),
                dominant_fraction=float(w.max() / w.sum()),
                transverse_mean_A=float(
                    np.average(projection["transverse_A"][mask], weights=w)
                ),
            )
        )
    return dict(
        histories=groups,
        support=support,
        direct_CV_identities_verified=3,
        exact_weight_order=True,
        sampling_converged=False,
        note="Empirical biased-history distributions; Kish ESS is not independent sampling ESS. No frames removed.",
    )
