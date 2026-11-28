"""Raw biased-lineage physical/proton checks, separate from PMF acceptance."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

NAMES = ("C_Nleave", "C_OG1", "OG1_HG1", "Nalpha_HG1", "Nleave_HG1")


def intersection(a: np.ndarray, b: np.ndarray) -> float:
    """Calculate histogram intersection with shared physical 0.05Å bins."""
    lo = np.floor(min(a.min(), b.min()) / 0.05) * 0.05
    hi = np.ceil(max(a.max(), b.max()) / 0.05) * 0.05 + 0.05
    edges = np.arange(lo, hi + 0.025, 0.05)
    x, y = np.histogram(a, edges)[0], np.histogram(b, edges)[0]
    return float(np.minimum(x / len(a), y / len(b)).sum())


def summarize(data: dict) -> dict:
    """Read exact retained frame identities for near/far and cross-start comparisons."""
    cache: dict[str, dict] = {}
    summaries = {}
    raw = {}
    offset = 0
    for state, count in zip(data["states"], data["counts"], strict=True):
        rows = []
        for sample in data["samples"][offset : offset + int(count)]:
            filename = sample["receipt"]
            if filename not in cache:
                receipt = Path(filename)
                proof = json.loads(receipt.read_text())
                with (receipt.parent / proof["csv_file"]).open() as stream:
                    cache[filename] = {
                        int(r["frame"]): r for r in csv.DictReader(stream)
                    }
            rows.append(cache[filename][sample["source_frame1based"]])
        values = np.asarray(
            [
                [
                    float(r["r_c_n"]),
                    float(r["r_c_og1"]),
                    float(r["r_og1_hg1"]),
                    float(r["r_og1_hg1"]) - float(r["q_act_A"]),
                    float(r["r_nleave_hg1"]),
                ]
                for r in rows
            ]
        )
        owner = [r["hg1_owner"] for r in rows]
        summary = {}
        for column, name in enumerate(NAMES):
            x = values[:, column]
            summary[name] = dict(
                mean=float(x.mean()),
                std=float(x.std()),
                first5_mean=float(x[:100].mean()),
                last5_mean=float(x[-100:].mean()),
                last_minus_first5=float(x[-100:].mean() - x[:100].mean()),
                first_last_histogram_intersection=intersection(x[:100], x[-100:]),
                five_ps_profiles=[
                    dict(
                        start_retained_ps=i * 0.05,
                        duration_ps=len(x[i : i + 100]) * 0.05,
                        mean=float(x[i : i + 100].mean()),
                    )
                    for i in range(0, len(x), 100)
                ],
            )
        summaries[state["state_id"]] = dict(
            distances=summary,
            hg1_owner_fractions={
                o: owner.count(o) / len(owner) for o in sorted(set(owner))
            },
            window_statistics=state["window_statistics"],
        )
        raw[state["state_id"]] = values
        offset += int(count)
    pairs = []
    for label, original in [
        ("A026_original", "activation_full_026"),
        ("A026_late", "activation_full_026"),
        ("A028_original", "activation_full_028"),
        ("A028_late", "activation_full_028"),
        ("D008_from_D010", "delivery_full_008"),
        ("D009_from_D007", "delivery_full_009"),
    ]:
        a, b = "main60/" + original, "mixing_checks30/" + label
        if a not in raw or b not in raw:
            raise ValueError(
                "Required independent-start lineage missing:" + a + ":" + b
            )
        pairs.append(
            dict(
                lineages=[a, b],
                histogram_intersections={
                    name: intersection(raw[a][:, i], raw[b][:, i])
                    for i, name in enumerate(NAMES)
                },
            )
        )
    return dict(
        lineages=summaries,
        independent_start_pairs=pairs,
        interpretation="Raw biased-state distributions; no pooling of independent starts. "
        "HG1 distances algebraically reconstructed from recorded coordinates. "
        "Nonstationarity flags remain; neither histogram overlap nor block CI proves convergence.",
    )
