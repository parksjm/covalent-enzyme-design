"""Fixed-weight conditional-bound sensitivity and dominant-frame receipt."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import logsumexp

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K"
)
SOURCE = BASE / "analysis/review_20261128_01/checkpoint_001_main70ps_supplement43"
OUT = Path(__file__).parent
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.relay_umbrella_cv import (  # noqa: E402
    RelayBiasWindow,
    amber_bias_energy,
    cross_bias_energies,
)
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


def write(name: str, value: Any) -> None:
    """Retain finite authoritative JSON evidence."""
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main() -> None:
    """Evaluate boundary and source-frame diagnostics without a weight refit."""
    data = np.load(SOURCE / "weights.npz")
    metadata = json.loads((SOURCE / "sources.json").read_text())
    q, lw, counts, fk = [data[k] for k in ["q", "logweights", "counts", "f_k"]]
    point = json.loads((SOURCE / "point.json").read_text())
    edges = np.arange(-2.5, 3.00001, 0.1)
    centers = (edges[1:] + edges[:-1]) / 2
    reference = int(
        np.searchsorted(edges, metadata["config"]["path_3cv_A"][0][1], side="right") - 1
    )
    index = 78151
    dominant_bin = int(np.searchsorted(edges, q[index, 1], side="right") - 1)
    outputs, rows = [], []
    for bound in [-2.0, -2.05, -2.10]:
        mask = (q[:, 0] >= bound) & (q[:, 0] < -0.8)
        h = PROJECTION.local_histogram(q[mask, 1], lw[mask], edges)
        ld = np.asarray(h["logdensity"], float)
        pmf = -KT * (ld - ld[reference])
        if bound == -2:
            np.testing.assert_allclose(
                pmf,
                np.asarray(point["projections"]["initial_q_act"]["pmf"], float),
                equal_nan=True,
            )
        local = (
            mask
            & (q[:, 1] >= edges[dominant_bin])
            & (q[:, 1] < edges[dominant_bin + 1])
        )
        fraction = (
            float(np.exp(lw[index] - logsumexp(lw[local]))) if mask[index] else 0.0
        )
        span = (centers >= -1.05 - 1e-9) & (centers <= 0.25 + 1e-9)
        use = np.flatnonzero(span)
        diffs = np.diff(pmf[use])
        outputs.append(
            dict(
                lower_bound_A=bound,
                upper_bound_A=-0.8,
                retained_frames=int(mask.sum()),
                dominant_frame_included=bool(mask[index]),
                dominant_qact_bin=float(centers[dominant_bin]),
                dominant_fraction_in_qact_bin=fraction,
                qact_bin_ess=float(h["ess"][dominant_bin]),
                monotonic_minus1p05_to_plus0p25=bool(np.all(diffs >= -1e-10)),
                maximum_downward_step_kcal_mol=float(max(0, -diffs.min())),
                pmf_at_minus0p85=float(pmf[dominant_bin]),
                pmf_at_plus0p25=float(pmf[use[-1]]),
                maximum_pmf_difference_from_original=float(
                    np.nanmax(
                        abs(
                            pmf[span]
                            - np.asarray(
                                point["projections"]["initial_q_act"]["pmf"], float
                            )[span]
                        )
                    )
                ),
            )
        )
        for i, x in enumerate(centers):
            if -1.5 <= x <= 0.8:
                local_bin = mask & (q[:, 1] >= edges[i]) & (q[:, 1] < edges[i + 1])
                rows.append(
                    dict(
                        lower_bound_A=bound,
                        q_act=float(x),
                        pmf_kcal_mol=float(pmf[i]),
                        raw_count=h["counts"][i],
                        kish_ess=h["ess"][i],
                        max_frame_fraction=float(
                            np.exp(lw[local_bin].max() - logsumexp(lw[local_bin]))
                        ),
                    )
                )
    write("sensitivity_summary.json", outputs)
    write("physical_profiles.json", rows)
    with (OUT / "physical_profiles.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    owner = np.repeat(np.arange(len(counts)), counts)
    state = metadata["states"][owner[index]]
    bias = state["bias"]
    deltas = q[index, :2] - np.asarray(bias["centers_A"])
    energies = [
        amber_bias_energy(float(v), c, k, bias["outer_extent_A"])
        for v, c, k in zip(
            q[index, :2], bias["centers_A"], bias["force_constants"], strict=True
        )
    ]
    windows = [
        RelayBiasWindow(
            tuple(x["bias"]["feature_names"]),
            tuple(x["bias"]["centers_A"]),
            tuple(x["bias"]["force_constants"]),
            x["bias"]["outer_extent_A"],
        )
        for x in metadata["states"]
    ]
    u = cross_bias_energies(
        dict(zip(["q_acyl", "q_act", "q_del"], q[index : index + 1].T, strict=True)),
        windows,
    )[:, 0]
    log_resp = np.log(counts) + fk - u / KT
    resp = np.exp(log_resp - logsumexp(log_resp))
    stage = (
        BASE / "activation_support_20261127/activation_support_3_03/04_production_5ps"
    )
    with (stage / "raw_features.csv").open() as handle:
        raw = list(csv.DictReader(handle))[51]
    np.testing.assert_allclose(
        [float(raw[k]) for k in ["q_acyl_A", "q_act_A", "q_del_A"]],
        q[index],
        rtol=0,
        atol=1e-14,
    )
    acceptance = json.loads((stage / "acceptance.json").read_text())
    analysis = json.loads((stage / "analysis_acceptance.json").read_text())
    receipt = dict(
        pooled_frame_index_zero_based=index,
        state_id=state["state_id"],
        q_A=q[index].tolist(),
        bias=bias,
        delta_from_bias_center_A=deltas.tolist(),
        bias_components_kcal_mol=energies,
        total_bias_kcal_mol=sum(energies),
        total_bias_in_kT=sum(energies) / KT,
        raw_feature_frame=int(raw["frame"]),
        raw_feature_row_zero_based=51,
        native_time_ps=float(raw["time_ps"]),
        production_age_ps=7.60,
        stage=str(stage),
        native_bias_kcal_mol=float(raw["native_bias_kcal_mol"]),
        native_bias_error_bound_kcal_mol=float(raw["native_bias_error_bound_kcal_mol"]),
        native_bias_evaluated=raw["native_bias_evaluated"],
        acceptance=acceptance,
        analysis_acceptance=analysis,
        source_state_qdel_range_A=[
            float(q[owner == owner[index], 2].min()),
            float(q[owner == owner[index], 2].max()),
        ],
        source_state_qdel_lower_rank=int(
            np.sum(q[owner == owner[index], 2] <= q[index, 2])
        ),
        nearest_bias_energy_kcal_mol=float(u.min()),
        nearest_bias_state=metadata["states"][int(u.argmin())]["state_id"],
        denominator_responsibilities=[
            dict(
                state=metadata["states"][i]["state_id"],
                fraction=float(resp[i]),
                bias_kcal_mol=float(u[i]),
            )
            for i in np.argsort(resp)[-5:][::-1]
        ],
    )
    write("dominant_frame_receipt.json", receipt)
    hard = project_path(q, np.asarray(metadata["config"]["path_3cv_A"]))
    w = np.exp(lw - logsumexp(lw))
    stats = []
    for left in [0.1, 0.15, 0.225, 0.25]:
        mask = (hard["s"] >= left) & (hard["s"] < left + 0.025)
        stats.append(
            dict(
                s=left + 0.0125,
                weighted_mean_transverse_A=float(
                    w[mask] @ hard["transverse_A"][mask] / w[mask].sum()
                ),
                minimum_transverse_A=float(hard["transverse_A"][mask].min()),
                dominant_transverse_A=float(
                    hard["transverse_A"][np.flatnonzero(mask)[np.argmax(w[mask])]]
                ),
            )
        )
    write("transverse_support.json", stats)
    inputs = [
        Path(__file__),
        SOURCE / "weights.npz",
        SOURCE / "sources.json",
        SOURCE / "point.json",
        HELPER,
        ROOT / "src/nyl12_covalent/relay_umbrella_cv.py",
        ROOT / "src/nyl12_covalent/relay_umbrella_projection.py",
        stage / "raw_features.csv",
        stage / "acceptance.json",
        stage / "analysis_acceptance.json",
    ]
    write(
        "input_sha256.json",
        {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    )
    write(
        "acceptance.json",
        dict(
            fixed_weights=True,
            baseline_reproduced=True,
            dominant_row_matches=True,
            sampling_converged=False,
        ),
    )
    logging.info("Boundary sensitivity and native dominant-frame receipt complete")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
