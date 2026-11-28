"""Strict whole-sequence audits and rankings for the bounded relay pilot."""

from __future__ import annotations

from typing import Any

import numpy as np

from .nyl12_sequence_protection import (
    MINIMUM_PROTECTED_ALPHA,
    MINIMUM_PROTECTED_BETA_LOCAL,
)

from .ligandmpnn_pilot_audit import (
    EXPECTED_CHAINS,
    _audit_ties_and_changes,
    _split_sequence,
)


def audit_sequence(
    sequence: str, native: dict[str, str], mask: dict[str, list[int]]
) -> dict[str, Any]:
    """Validate one complete tied sequence and return canonical mutations."""
    chains = _split_sequence(sequence, "relay_ensemble", 0)
    changes, alpha, beta = _audit_ties_and_changes(
        chains,
        native,
        set(mask["alpha"]),
        set(mask["beta"]),
        set(MINIMUM_PROTECTED_ALPHA),
        set(MINIMUM_PROTECTED_BETA_LOCAL),
        "relay_ensemble",
        0,
    )
    return {
        "sequence": sequence,
        "alpha_sequence": alpha,
        "beta_sequence": beta,
        "mutations": changes,
        "mutation_count": len(changes),
    }


def pool_samples(samples: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Deduplicate entire sequences, retaining every sampling lineage."""
    pool: dict[str, dict[str, Any]] = {}
    for sample in samples:
        key = sample["sequence"]
        if key not in pool:
            pool[key] = {
                **sample,
                "candidate_id": f"relay_{len(pool) + 1:02d}",
                "lineage": [],
            }
        pool[key]["lineage"].append(sample["sample_id"])
    if len(pool) > 24:
        raise ValueError("Unique pool exceeds the authorized 24-sequence budget")
    return list(pool.values())


def summarize_score(
    stats: dict[str, Any], sequence: str, mask: dict[str, list[int]]
) -> dict[str, float]:
    """Gather candidate log likelihoods with whole, active and per-copy means."""
    chains = _split_sequence(sequence, "crossscore", 0)
    names = stats["residue_names"]
    names = [names[i] for i in range(len(names))]
    expected = [
        f"{chain}{p}" for chain, n in EXPECTED_CHAINS.items() for p in range(1, n + 1)
    ]
    if names != expected:
        raise ValueError("Native score residue order differs from canonical mapping")
    alphabet = stats["alphabet"]
    indices = np.array(
        [alphabet.index(aa) for c in EXPECTED_CHAINS for aa in chains[c]]
    )
    if not np.array_equal(indices, stats["native_sequence"]):
        raise ValueError("Native scoring did not use the complete candidate sequence")
    log_probs = np.asarray(stats["log_probs"])
    if log_probs.ndim != 3 or log_probs.shape[1:] != (1312, len(alphabet)):
        raise ValueError("Unexpected native score dimensions")
    if not np.isfinite(log_probs).all() or not np.all(stats["mask"]):
        raise ValueError("Nonfinite scores or missing backbone residues")
    nll = -log_probs[:, np.arange(1312), indices].mean(axis=0)
    result = {"whole_nll": float(nll.mean())}
    all_designed = []
    for site, pair in zip(("AE", "BF", "CG", "DH"), zip("ABCD", "EFGH")):
        chosen = [
            i
            for i, name in enumerate(names)
            if name[0] in pair
            and int(name[1:]) in mask["alpha" if name[0] in "ABCD" else "beta"]
        ]
        if not chosen:
            raise ValueError("Empty designed-position score mask")
        result[f"{site}_design_nll"] = float(nll[chosen].mean())
        all_designed.extend(chosen)
    result["allcopy_design_nll"] = float(nll[all_designed].mean())
    return result


def rank_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rank by active-template mean and flag the mean/worst Pareto frontier."""
    for candidate in candidates:
        values = [s["AE_design_nll"] for s in candidate["scores"] if s["state"] == "TS"]
        if len(values) != 3 or not np.isfinite(values).all():
            raise ValueError("Every candidate requires three finite TS crossscores")
        candidate["ts_mean_nll"] = float(np.mean(values))
        candidate["ts_worst_nll"] = float(np.max(values))
    ordered = sorted(candidates, key=lambda c: (c["ts_mean_nll"], c["candidate_id"]))
    for rank, candidate in enumerate(ordered, 1):
        candidate["rank"] = rank
        candidate["pareto"] = not any(
            other["ts_mean_nll"] <= candidate["ts_mean_nll"]
            and other["ts_worst_nll"] <= candidate["ts_worst_nll"]
            and other["mutation_count"] <= candidate["mutation_count"]
            and (
                other["ts_mean_nll"] < candidate["ts_mean_nll"]
                or other["ts_worst_nll"] < candidate["ts_worst_nll"]
                or other["mutation_count"] < candidate["mutation_count"]
            )
            for other in ordered
        )
    return ordered
