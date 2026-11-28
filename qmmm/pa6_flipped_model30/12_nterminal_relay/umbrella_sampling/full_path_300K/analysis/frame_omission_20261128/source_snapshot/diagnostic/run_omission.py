"""Recompute a one-frame-omitted PMF for curiosity; never production analysis."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import logging
from pathlib import Path
import sys
import time

import numpy as np
from scipy.special import logsumexp

from omission_core import exact_refit, omit_observation, profile
from omission_figures import plot_profiles
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from nyl12_covalent.relay_umbrella_projection import project_path

KT = 0.00198720425864083 * 300
LOGGER = logging.getLogger(__name__)


def encode(value: object) -> object:
    """Convert numpy arrays and unsupported bins to finite JSON values."""
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def write(path: Path, value: dict) -> None:
    """Write authoritative finite JSON evidence."""
    path.write_text(json.dumps(encode(value), indent=2, allow_nan=False) + "\n")


def compare(result: dict, baseline: dict, centers: np.ndarray, target: int) -> dict:
    """Summarize local concentration and global drift on common supported bins."""
    common = np.isfinite(result["pmf"]) & np.isfinite(baseline["pmf"])
    delta = result["pmf"][common] - baseline["pmf"][common]
    outside = common.copy()
    outside[target] = False
    peak = int(np.nanargmax(result["pmf"]))
    base_peak = int(np.nanargmax(baseline["pmf"]))
    return dict(
        affected_bin_s=float(centers[target]),
        affected_bin_pmf_kcal_mol=float(result["pmf"][target]),
        affected_bin_change_kcal_mol=float(
            result["pmf"][target] - baseline["pmf"][target]
        ),
        affected_bin_count=int(result["counts"][target]),
        affected_bin_kish_ess=float(result["ess"][target]),
        affected_bin_max_fraction=float(result["max_fraction"][target]),
        affected_bin_dominant_original_index=int(
            result["dominant_original_index"][target]
        ),
        maximum_occupied_bin_kcal_mol=float(result["pmf"][peak]),
        maximum_occupied_bin_s=float(centers[peak]),
        peak_height_change_kcal_mol=float(
            result["pmf"][peak] - baseline["pmf"][base_peak]
        ),
        common_supported_bins=int(common.sum()),
        lost_supported_bins=np.flatnonzero(
            np.isfinite(baseline["pmf"]) & ~common
        ).tolist(),
        common_bin_max_abs_change_kcal_mol=float(abs(delta).max()),
        common_bin_rms_change_kcal_mol=float(np.sqrt(np.mean(delta**2))),
        common_bin_median_change_kcal_mol=float(np.median(delta)),
        common_bin_min_change_kcal_mol=float(delta.min()),
        common_bin_max_change_kcal_mol=float(delta.max()),
        outside_affected_bin_max_abs_change_kcal_mol=float(
            np.max(abs(result["pmf"][outside] - baseline["pmf"][outside]))
        ),
    )


def main() -> None:
    """Verify frozen provenance, refit one omission and retain a diagnostic report."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    base = Path(sys.argv[1])
    inputs, output = base / "frozen_inputs", base / "attempt_01/results"
    output.mkdir(exist_ok=False, parents=True)
    manifest = json.loads((base / "provenance.json").read_text())
    for name, expected in manifest["sha256"].items():
        if hashlib.sha256((base / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Frozen source/input hash mismatch: {name}")
    data = np.load(inputs / "weights.npz")
    q, counts, lw, fk = [data[k] for k in ("q", "counts", "logweights", "f_k")]
    sources = json.loads((inputs / "sources.json").read_text())
    point = json.loads((inputs / "point.json").read_text())
    receipt = json.loads((inputs / "dominant_frame_receipt.json").read_text())
    if q.shape != (82000, 3) or counts.shape != (100,) or counts.sum() != 82000:
        raise ValueError(
            "Expected exactly the immutable 82000-frame, 100-history baseline"
        )
    if sources["config"]["temperature_K"] != 300:
        raise ValueError("Expected original 300 K analysis")
    states = sources["states"]
    state_ids = [
        i for i, st in enumerate(states) if st["state_id"] == receipt["state_id"]
    ]
    if len(state_ids) != 1:
        raise ValueError("Source state ID must be unique")
    state = state_ids[0]
    with (inputs / "raw_features.csv").open() as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if int(row["frame"]) == receipt["raw_feature_frame"]
            and abs(float(row["time_ps"]) - receipt["native_time_ps"]) < 1e-10
        ]
    if len(matches) != 1:
        raise ValueError("Native stage/frame/time must identify exactly one source row")
    target_q = np.array(
        [float(matches[0][f"{name}_A"]) for name in ("q_acyl", "q_act", "q_del")]
    )
    np.testing.assert_allclose(target_q, receipt["q_A"], atol=1e-14, rtol=0)
    expected_csv = receipt["analysis_acceptance"]["csv_sha256"]
    if (
        hashlib.sha256((inputs / "raw_features.csv").read_bytes()).hexdigest()
        != expected_csv
    ):
        raise ValueError("Native source does not match the accepted CSV digest")
    keep, new_counts, index = omit_observation(q, counts, state, target_q)
    if index != receipt["pooled_frame_index_zero_based"]:
        raise ValueError(
            "Provenance-and-coordinate match disagrees with historical pooled ID"
        )
    local_index = index - int(counts[:state].sum())
    selected = states[state]["window_statistics"]["selected_indices"][local_index]
    production_age = (selected + 1) * states[state]["window_statistics"][
        "frame_interval_ps"
    ]
    np.testing.assert_allclose(production_age, receipt["production_age_ps"], atol=1e-12)
    write(
        output / "omitted_observation.json",
        dict(
            diagnostic_only=True,
            original_receipt=receipt,
            matched_pooled_index=index,
            owner_state_index=state,
            owner_local_retained_index=local_index,
            selected_history_row=selected,
            production_age_ps=production_age,
            original_N_k=int(counts[state]),
            omitted_N_k=int(new_counts[state]),
            observation_removed_from_analysis_only=True,
            original_files_preserved=True,
        ),
    )
    config = sources["config"]
    edges = np.array(config["s_bin_edges"])
    centers = (edges[:-1] + edges[1:]) / 2
    np.testing.assert_allclose(np.diff(edges), 0.025, atol=1e-12)
    reference = int(np.argmin(abs(centers - 0.0125)))
    np.testing.assert_allclose(centers[reference], 0.0125, atol=1e-12)
    s = project_path(q, np.array(config["path_3cv_A"]))["s"]
    target = int(np.searchsorted(edges, s[index], side="right") - 1)
    original_indices = np.arange(len(q))
    baseline = profile(s, lw, edges, reference, KT, original_indices)
    parity = float(
        np.nanmax(
            abs(baseline["pmf"] - np.array(point["projections"]["hard"]["pmf"], float))
        )
    )
    if parity > 1e-8:
        raise ValueError(f"Baseline point parity failure: {parity}")
    windows = [RelayBiasWindow(**st["bias"]) for st in states]
    u = (
        cross_bias_energies(
            dict(zip(("q_acyl", "q_act", "q_del"), q.T, strict=True)), windows
        )
        / KT
    )
    reconstructed = -logsumexp(np.log(counts)[:, None] + fk[:, None] - u, axis=0)
    reconstructed -= logsumexp(reconstructed)
    weight_parity = float(abs(reconstructed - lw).max())
    if weight_parity > 1e-7:
        raise ValueError(f"Native cross-bias weight parity failure: {weight_parity}")
    LOGGER.info(
        "Verified single observation %s in state %s; starting 81999-frame MBAR",
        index,
        state,
    )
    start = time.monotonic()
    fit = exact_refit(u[:, keep], new_counts, fk)
    full = profile(
        s[keep], fit["logweights"], edges, reference, KT, original_indices[keep]
    )
    fixed = profile(s[keep], lw[keep], edges, reference, KT, original_indices[keep])
    share = float(
        np.exp(
            lw[index] - logsumexp(lw[(s >= edges[target]) & (s < edges[target + 1])])
        )
    )
    analytical = float(-KT * np.log1p(-share))
    np.testing.assert_allclose(
        fixed["pmf"][target] - baseline["pmf"][target], analytical, atol=1e-10
    )
    curves = dict(baseline=baseline, full_refit=full, fixed_weights=fixed)
    for curve in curves.values():
        np.testing.assert_allclose(curve["pmf"][reference], 0, atol=1e-12)
        np.testing.assert_allclose(curve["mass"].sum(), 1, atol=1e-12)
    summary = {
        key: compare(value, baseline, centers, target) for key, value in curves.items()
    }
    new_dom = summary["full_refit"]["affected_bin_dominant_original_index"]
    owner = np.repeat(np.arange(len(counts)), counts)
    result = dict(
        label="DIAGNOSTIC ONLY — NOT FOR PUBLICATION",
        sampling_converged=False,
        production_primary_unchanged=True,
        no_physical_validity_rejection_of_frame=True,
        interpretation="Single-observation sensitivity only; no accepted activation barrier or statistical confidence interval",
        frames_before=len(q),
        frames_after=int(keep.sum()),
        states=len(counts),
        counts_before=counts,
        counts_after=new_counts,
        reference_bin=reference,
        reference_s=centers[reference],
        bin_width_s=0.025,
        temperature_K=300,
        kt_kcal_mol=KT,
        native_cross_bias_logweight_max_error=weight_parity,
        baseline_point_max_error_kcal_mol=parity,
        fixed_weight_analytical_shift_kcal_mol=analytical,
        original_affected_bin_dominant_share=share,
        summary=summary,
        curves=curves,
        s_centers=centers,
        s_edges=edges,
        new_dominant_frame=dict(
            original_pooled_index=new_dom,
            state_id=states[owner[new_dom]]["state_id"],
            q_A=q[new_dom],
        ),
        fit={
            k: v
            for k, v in fit.items()
            if k not in ("logweights", "overlap_matrix", "f_k")
        },
        elapsed_refit_seconds=time.monotonic() - start,
        versions={
            key: importlib.metadata.version(key)
            for key in ("numpy", "scipy", "pymbar", "matplotlib")
        },
        python=sys.executable,
    )
    write(output / "comparison.json", result)
    np.savez_compressed(
        output / "refitted_weights.npz",
        original_indices=original_indices[keep],
        q=q[keep],
        counts=new_counts,
        f_k=fit["f_k"],
        logweights=fit["logweights"],
        fixed_original_logweights=lw[keep] - logsumexp(lw[keep]),
        overlap_matrix=fit["overlap_matrix"],
    )
    plot_profiles(centers, curves, output)
    rows = []
    for i, center in enumerate(centers):
        rows.append(
            dict(
                s=center,
                **{
                    f"{key}_{field}": value[field][i]
                    for key, value in curves.items()
                    for field in ("pmf", "ess", "max_fraction", "counts")
                },
            )
        )
    with (output / "profiles.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write(
        output / "scientific_acceptance.json",
        dict(
            diagnostic_checks_passed=True,
            sampling_converged=False,
            accepted_activation_barrier=False,
            primary_analysis_unchanged=True,
            process_completion="See ../status.json",
            checks=[
                "Accepted native CSV digest plus state/frame/time/CV match",
                "Exactly one observation omitted and only its N_k decremented",
                "Frozen 0.025 bins and s=0.0125 zero",
                "Baseline PMF and native cross-bias logweight parity",
                "Robust full 100-state MBAR refit with gauge residual <=1e-7 kBT",
                "Overlap graph connected at 1e-12",
                "Fixed-weight local shift matches -kT log(1-share)",
                "All original sampling files preserved",
            ],
        ),
    )
    report = [
        "# DIAGNOSTIC ONLY — NOT FOR PUBLICATION",
        "",
        "Curiosity-only omission of one accepted observation; production MBAR is unchanged.",
        "Full refit: all native bias terms, 81,999 observations, 100 states; decrement only state 85 from 100 to 99.",
        "Fixed original weights: delete the same observation and renormalize; no MBAR refit.",
        "Both use the original 0.025-wide s bins and zero at s=0.0125. No new uncertainty or accepted barrier.",
        "",
        "| Estimate | F at s=0.2625 | ESS there | Maximum frame share | Whole-profile maximum |",
        "|---|---:|---:|---:|---:|",
    ]
    for key, val in summary.items():
        report.append(
            f"| {key} | {val['affected_bin_pmf_kcal_mol']:.9f} | {val['affected_bin_kish_ess']:.6f} | "
            f"{100 * val['affected_bin_max_fraction']:.5f}% | {val['maximum_occupied_bin_kcal_mol']:.9f} at s={val['maximum_occupied_bin_s']:.4f} |"
        )
    report += [
        "",
        "Energies are kcal/mol. Full authoritative evidence: comparison.json and omitted_observation.json.",
        f"Fixed-weight analytical local shift: {analytical:.9f} kcal/mol.",
        f"Full-refit maximum common-bin change: {summary['full_refit']['common_bin_max_abs_change_kcal_mol']:.9f} kcal/mol; "
        f"outside affected bin: {summary['full_refit']['outside_affected_bin_max_abs_change_kcal_mol']:.9f} kcal/mol.",
        f"Gauge residual: {fit['fixed_point_residual_kBT']:.3g} kBT; overlap components: {fit['connected_components_at_1e_minus12']}.",
        f"New dominant observation in affected bin: {new_dom} ({states[owner[new_dom]]['state_id']}).",
        "",
        "This diagnoses sensitivity to one high-weight accepted observation. Its removal is not justified as data cleaning.",
        "Neither estimator establishes sampling convergence, a physical transition state or a validated activation barrier.",
    ]
    (output / "README.md").write_text("\n".join(report) + "\n")
    LOGGER.info("Completed diagnostic: %s", summary)


if __name__ == "__main__":
    main()
