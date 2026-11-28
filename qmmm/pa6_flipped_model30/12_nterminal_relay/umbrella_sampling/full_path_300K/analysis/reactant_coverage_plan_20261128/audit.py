"""Audit exact MBAR coverage without reading raw trajectories or changing jobs."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from nyl12_covalent.relay_umbrella_projection import project_path

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K"
)
OUT = Path(__file__).resolve().parent
SOURCE = (
    BASE
    / "analysis/main100_continuation_20261127/current_ext100/checkpoint_032_main70ps_supplement43"
)


def main() -> None:
    """Write source-hashed per-state responsibility and neighborhood diagnostics."""
    source = json.loads((SOURCE / "sources.json").read_text())
    data = np.load(SOURCE / "weights.npz")
    q, counts, free = data["q"], data["counts"], data["f_k"]
    assert q.shape == (82000, 3) and len(counts) == 100
    states = source["states"]
    feature_ids = {"q_acyl": 0, "q_act": 1, "q_del": 2}
    beta = 1 / (0.00198720425864083 * source["config"]["temperature_K"])
    energies = np.zeros((len(states), len(q)))
    windows = []
    for i, state in enumerate(states):
        bias = state["bias"]
        windows.append(RelayBiasWindow(**bias))
        for name, center, force in zip(
            bias["feature_names"],
            bias["centers_A"],
            bias["force_constants"],
            strict=True,
        ):
            delta = np.abs(q[:, feature_ids[name]] - center)
            extent = bias["outer_extent_A"]
            energies[i] += force * np.where(
                delta <= extent, delta**2, extent * (2 * delta - extent)
            )
    check = cross_bias_energies(
        {name: q[:10, col] for name, col in feature_ids.items()}, windows
    )
    assert np.allclose(energies[:, :10], check, atol=1e-12)
    terms = np.log(counts[:, None]) + free[:, None] - beta * energies
    denominator = logsumexp(terms, axis=0)
    normalized = -denominator - logsumexp(-denominator)
    assert np.max(np.abs(normalized - data["logweights"])) < 1e-5
    responsibilities = np.exp(terms - denominator)
    point = q[78151]
    projection = project_path(q, np.array(source["config"]["path_3cv_A"]))
    s = projection["s"]
    local3 = np.all(np.abs(q - point) <= 0.15, axis=1)
    local2 = np.all(np.abs(q[:, :2] - point[:2]) <= 0.15, axis=1)
    region = (q[:, 0] < -1.9) & (q[:, 1] > -1.1) & (q[:, 1] < -0.5) & (q[:, 2] < -2.0)
    sbin = (s >= 0.25) & (s < 0.275)
    offsets = np.r_[0, np.cumsum(counts)]
    rows = []
    hours = []
    live = []
    for i, state in enumerate(states):
        block = q[offsets[i] : offsets[i + 1]]
        masks = {
            "within_0p15_A_3CV": local3,
            "within_0p15_A_2CV": local2,
            "reactant_tail_region": region,
            "s_0p25_0p275": sbin,
        }
        row = {
            "state_id": state["state_id"],
            "counts": int(counts[i]),
            "bias": state["bias"],
            "responsibility_at_dominant": float(responsibilities[i, 78151]),
            "bias_at_dominant_kcal_mol": float(energies[i, 78151]),
            "q_mean": block.mean(axis=0).tolist(),
            "q_std": block.std(axis=0).tolist(),
            "q_quantiles_01_10_50_90_99": np.quantile(
                block, [0.01, 0.1, 0.5, 0.9, 0.99], axis=0
            ).tolist(),
            "neighborhood_counts": {
                name: int(mask[offsets[i] : offsets[i + 1]].sum())
                for name, mask in masks.items()
            },
            "g_by_column_diagnostic": state["window_statistics"]["g_by_column"],
        }
        rows.append(row)
        if state["group"] == "activation_grid33":
            folder = BASE / "activation_support_20261127" / state["window"]
            processes = []
            for file in sorted(folder.glob("*/process.json")):
                record = json.loads(file.read_text())
                if "ended_utc" in record and "production" in file.parent.name:
                    elapsed = (
                        datetime.fromisoformat(record["ended_utc"])
                        - datetime.fromisoformat(record["started_utc"])
                    ).total_seconds() / 3600
                    hours.append(elapsed)
                processes.append(
                    {
                        "stage": file.parent.name,
                        **{
                            key: record[key]
                            for key in ["started_utc", "ended_utc", "exit_status"]
                            if key in record
                        },
                    }
                )
            live.append(
                {
                    "state_id": state["state_id"],
                    "accepted_production_horizon_ps": state["production_horizon_ps"],
                    "processes": processes,
                }
            )
    totals = {name: int(mask.sum()) for name, mask in masks.items()}
    ordered = sorted(rows, key=lambda row: -row["responsibility_at_dominant"])
    group_resp = {
        group: float(
            sum(
                row["responsibility_at_dominant"]
                for row in rows
                if row["state_id"].startswith(group + "/")
            )
        )
        for group in {state["group"] for state in states}
    }
    result = {
        "source": str(SOURCE),
        "source_sha256": {
            name: hashlib.sha256((SOURCE / name).read_bytes()).hexdigest()
            for name in ["weights.npz", "sources.json"]
        },
        "dominant_index": 78151,
        "dominant_q_A": point.tolist(),
        "bias_convention": "Amber k*(q-center)^2, tangent-linear beyond extent; responsibilities=N_k exp(f_k-beta V_kn)/sum_j",
        "same_method_check_max_abs_logweight_error": float(
            np.max(np.abs(normalized - data["logweights"]))
        ),
        "neighborhood_totals": totals,
        "responsibility_by_group": group_resp,
        "top_responsibilities": ordered[:10],
        "all_states": rows,
        "sigma_A_for_k25": float(np.sqrt(1 / (beta * 50))),
        "timing_completed_production_5ps_h": {
            "n": len(hours),
            "quantiles_10_50_90": np.quantile(hours, [0.1, 0.5, 0.9]).tolist(),
        },
        "live_support_at_audit": live,
    }
    (OUT / "coverage_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    summary = {
        key: result[key]
        for key in [
            "dominant_q_A",
            "neighborhood_totals",
            "responsibility_by_group",
            "timing_completed_production_5ps_h",
        ]
    }
    summary["top_responsibilities"] = [
        {
            key: row[key]
            for key in [
                "state_id",
                "responsibility_at_dominant",
                "bias_at_dominant_kcal_mol",
                "neighborhood_counts",
            ]
        }
        for row in ordered[:6]
    ]
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
