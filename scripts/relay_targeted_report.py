"""Build a bounded recommendation packet from accepted targeted trajectories."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from relay_targeted_statistics import (
    FEATURES,
    histogram_overlap,
    production_rows,
    recommendation,
    require_same_bias,
    summarize,
)


def make_report(base: Path, cfg: dict, output: Path) -> dict[str, Any]:
    """Write JSON, Markdown and unsmoothed relaxation/histogram evidence.

    Args:
        base: Existing parent full-path campaign directory.
        cfg: Frozen four-target preparation.
        output: Existing targeted analysis directory.

    Returns:
        Authoritative diagnostic report.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parent = json.loads((base / "prepared.json").read_text())
    windows = {r["label"]: r for r in parent["windows"]}
    data, report = {}, {"windows": {}, "comparisons": {}}
    for row in cfg["windows"]:
        label = row["label"]
        paths = [
            base / "additional_sampling_20261126" / s["input"]
            for s in row["stages"]
            if s["kind"] == "production"
        ]
        data[label] = production_rows([p.parent / "raw_features.csv" for p in paths])
        report["windows"][label] = dict(
            summary=summarize(data[label][1]),
            lineage=row["source_lineage"],
            fresh_velocity_input=(
                base / "additional_sampling_20261126" / row["stages"][0]["input"]
            ).read_text(),
            seeds=[s["seed"] for s in row["stages"]],
            bias=row["bias"],
        )
    for label in ("activation_full_021", "activation_full_022", "delivery_full_026"):
        paths = [
            base / label / f"{i:02}_production_5ps" / "raw_features.csv"
            for i in range(3, 9)
        ]
        data[label] = production_rows(paths)
        report["windows"][label] = dict(
            summary=summarize(data[label][1]), bias=windows[label]["bias"]
        )
    pairs = [
        (r["label"], r["source_lineage"]["opposite_start_target"], True)
        for r in cfg["windows"]
        if r["label"].startswith("activation")
    ]
    products = ["delivery_full_026", "delivery_027_qdel_1p50", "delivery_028_qdel_1p75"]
    pairs += [(a, b, False) for a, b in zip(products, products[1:])]
    for a, b, same in pairs:
        if same:
            require_same_bias(
                report["windows"][a]["bias"], report["windows"][b]["bias"]
            )

        def compare(ra: list[dict], rb: list[dict]) -> dict:
            values = {}
            for feature in FEATURES:
                values[feature] = histogram_overlap(
                    np.array([float(r[feature]) for r in ra]),
                    np.array([float(r[feature]) for r in rb]),
                )
            cvs = ["q_acyl_A", "q_act_A" if same else "q_del_A"]
            values["joint_biased_CV"] = histogram_overlap(
                np.array([[float(r[k]) for k in cvs] for r in ra]),
                np.array([[float(r[k]) for k in cvs] for r in rb]),
            )
            return values

        overlap = compare(data[a][1], data[b][1])
        stats = [report["windows"][x]["summary"]["features"] for x in (a, b)]
        ess = [v["ess"] for s in stats for v in s.values()]
        gap = abs(stats[0]["q_del_A"]["mean"] - stats[1]["q_del_A"]["mean"])
        drift = max(abs(v["half_drift_A"]) for s in stats for v in s.values())
        report["comparisons"][a + "__" + b] = dict(
            same_bias=same,
            histogram_intersection=overlap,
            matched_last5ps=compare(
                *[
                    [
                        r
                        for r in data[x][0]
                        if r["production_time_ps"]
                        > float(data[x][0][-1]["production_time_ps"]) - 5 + 1e-8
                    ]
                    for x in (a, b)
                ]
            ),
            qdel_mean_gap_A=gap,
            max_half_drift_A=drift,
            recommendation=recommendation(
                min(ess) if all(v is not None for v in ess) else None,
                drift,
                overlap["q_del_A"],
                same,
                overlap["joint_biased_CV"],
                gap,
            ),
        )
    report["protocol"] = dict(
        burnin="Discard first 5 production ps once per independent trajectory; exclude all equilibration",
        histogram="Pooled zero-anchored 0.1 A bins; no censoring or smoothing",
        thresholds=dict(
            low_overlap=0.2,
            same_bias_agreement=0.5,
            qdel_gap_A=0.3,
            ess_warning=20,
            half_drift_A=0.2,
        ),
        limitations="Prospective diagnostic heuristics, not universal proof. Short post-burnin samples and ESS uncertainty limit decisions. No PMF or convergence claim. Product tests need >=30 ps before final matched PMF claims. Chemistry retention requires lead review of distances and full proton owners; no automatic chemistry classification.",
        next_action="Recommendation only; never launch windows/extensions automatically",
    )
    report["chemistry_acceptance"] = "requires_lead_review"
    report["all_recommendations_conditional"] = True
    friendly = dict(
        zip(
            data,
            [
                "21 ← 22",
                "22 ← 21",
                "27 diagnostic",
                "28 diagnostic",
                "21 reference",
                "22 reference",
                "26 reference",
            ],
        )
    )
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 14, "legend.fontsize": 10})
    for group in ("activation", "delivery"):
        fig, axes = plt.subplots(2, 2, figsize=(13, 8))
        for label, (full, kept) in data.items():
            if not label.startswith(group):
                continue
            times = [r["production_time_ps"] for r in full]
            axes[0, 0].plot(
                times,
                [float(r["q_del_A"]) for r in full],
                lw=0.7,
                label=friendly[label],
            )
            x = np.array([float(r["q_del_A"]) for r in kept])
            bins = np.arange(
                np.floor(x.min() * 10) / 10, np.ceil(x.max() * 10) / 10 + 0.2, 0.1
            )
            axes[0, 1].hist(
                x, bins=bins, histtype="step", density=True, label=friendly[label]
            )
            for ax, key, title in zip(
                axes[1],
                ("r_c_og1", "r_c_n"),
                ("C–Oγ1 distance (Å)", "C–N distance (Å)"),
            ):
                ax.plot(times, [float(r[key]) for r in full], lw=0.7)
                ax.set_ylabel(title)
                ax.set_xlabel("Production time (ps)")
        axes[0, 0].set_ylabel("$q_{del}$ (Å)")
        axes[0, 0].set_xlabel("Production time (ps)")
        axes[0, 1].set_xlabel("$q_{del}$ after fixed burn-in (Å)")
        axes[0, 1].set_ylabel("Probability density (Å⁻¹)")
        axes[0, 0].legend()
        fig.suptitle(group.capitalize() + " diagnostic sampling")
        fig.tight_layout()
        fig.savefig(output / f"targeted_relay_{group}_relaxation.png", dpi=160)
        plt.close(fig)
    (output / "targeted_relay_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    lines = [
        "# Targeted relay diagnostic screens",
        "",
        report["protocol"]["limitations"],
        "",
        "| Comparison | qdel overlap | gap A | Recommendation |",
        "|---|---:|---:|---|",
    ]
    for pair, c in report["comparisons"].items():
        lines.append(
            f"| {pair} | {c['histogram_intersection']['q_del_A']:.3f} | {c['qdel_mean_gap_A']:.3f} | {c['recommendation']} |"
        )
    lines += [
        "",
        report["protocol"]["burnin"],
        "",
        "See JSON for every raw-frame statistic, proton inventory, ESS, late-half drift, lineage, seeds and matched-last-5-ps sensitivity. No frames were censored. Inspect chemistry before accepting further sampling.",
    ]
    (output / "targeted_relay_report.md").write_text("\n".join(lines) + "\n")
    return report
