"""Tests for shared fixed-window umbrella statistics."""

from __future__ import annotations

import numpy as np

from nyl12_covalent.fixed_us_statistics import (
    _features,
    fit_fixed_data,
    run_fixed_statistics,
)
from nyl12_covalent.neb_umbrella_analysis import KB
from nyl12_covalent.neb_umbrella_bias import HarmonicWindow, bias_energy
from nyl12_covalent.relay_fixed_burnin import fixed_window_data
from nyl12_covalent.relay_fixed_pmf import fit_samples
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow, cross_bias_energies
from nyl12_covalent.relay_umbrella_projection import project_path


def _bias(center: float) -> dict:
    return dict(
        feature_names=("q_acyl", "q_act"),
        centers_A=(center, 0.0),
        force_constants=(1.0, 1.0),
    )


def _fixture() -> tuple[dict, list[dict]]:
    rng = np.random.default_rng(17)
    config = dict(
        temperature_K=300.0,
        windows=[dict(label="a", bias=_bias(-0.2)), dict(label="b", bias=_bias(0.2))],
        path_3cv_A=[[-0.6, 0.0, 0.0], [0.0, 0.0, 0.0], [0.6, 0.0, 0.0]],
        s_bin_edges=np.linspace(-0.5, 1.5, 9),
        fixed_reference_s=0.125,
        sample_feature_names=("q_acyl", "q_act", "q_del"),
    )
    blocks = []
    for label, center in (("a", -0.2), ("b", 0.2)):
        values = np.column_stack(
            (
                rng.normal(center, 0.18, 40),
                rng.normal(0.0, 0.08, 40),
                rng.normal(0.0, 0.08, 40),
            )
        )
        times = 0.5 * np.arange(1, 41)
        blocks.append(
            dict(
                window=label,
                values=values,
                original_rows=[dict(time_ps=str(time)) for time in times],
                sampled_duration_ps=20.0,
            )
        )
    return config, blocks


def test_derived_qh_keeps_native_direct_bias_energy_and_names() -> None:
    values = np.array([[0.2, 0.3, -0.1], [-0.4, -0.2, 0.05]])
    config = dict(
        sample_feature_names=("q_acyl", "q_act", "q_del"),
        windows=[dict(bias=dict(feature_names=("q_acyl", "q_h")))],
    )
    features = _features(config, values)
    np.testing.assert_allclose(features["q_h"], values[:, 1] + values[:, 2])
    window = RelayBiasWindow(
        feature_names=("q_acyl", "q_h"),
        centers_A=(0.1, -0.2),
        force_constants=(100.0, 100.0),
    )
    actual = cross_bias_energies(features, [window])[0]
    legacy = np.array(
        [
            bias_energy(
                (row[0], row[1] + row[2]),
                HarmonicWindow(0.1, -0.2, 100.0, 100.0),
            )
            for row in values
        ]
    )
    np.testing.assert_allclose(actual, legacy, rtol=0, atol=1e-12)
    assert "q_h" not in config["sample_feature_names"]


def test_fit_fixed_data_matches_existing_relay_point_estimate() -> None:
    config, blocks = _fixture()
    actual, _, summaries = fit_fixed_data(config, blocks, 5.0, 15.0)
    values, counts, expected_summaries = fixed_window_data(
        config["windows"], blocks, burn_in_ps=5.0, end_ps=15.0
    )
    biases = [RelayBiasWindow(**window["bias"]) for window in config["windows"]]
    features = dict(zip(("q_acyl", "q_act", "q_del"), values[:, :3].T, strict=True))
    reduced = cross_bias_energies(features, biases) / (KB * config["temperature_K"])
    projected = project_path(values[:, :3], np.asarray(config["path_3cv_A"]))["s"]
    expected = fit_samples(
        reduced,
        counts,
        projected,
        np.asarray(config["s_bin_edges"]),
        config["temperature_K"],
        config["fixed_reference_s"],
    )
    expected.update(burn_in_ps=5.0, end_ps=15.0)
    assert actual == expected
    assert summaries == expected_summaries


def test_fixed_statistics_preserves_bootstrap_support_semantics() -> None:
    config, blocks = _fixture()
    result = run_fixed_statistics(config, blocks, 15.0, replicates=50, workers=1)
    assert result["primary"]["connected"]
    assert set(result["bootstrap"]) == {"5ps", "2.5ps"}
    for bootstrap in result["bootstrap"].values():
        assert bootstrap["replicates"] == 50
        for support, lower, upper in zip(
            bootstrap["support_fraction"],
            bootstrap["lower_95_kcal_mol"],
            bootstrap["upper_95_kcal_mol"],
            strict=True,
        ):
            if support < 1.0:
                assert lower is None and upper is None
            else:
                assert lower is not None and upper is not None
    assert set(result["sensitivity"]) == {
        "burn_0ps",
        "burn_2.5ps",
        "burn_10ps",
        "early_half",
        "late_half",
    }
    assert [row["selected_frames"] for row in result["window_selection"]] == [
        20,
        20,
    ]


def test_empty_outer_bins_leave_common_pmf_bins_unchanged() -> None:
    """Grid coverage can expand without changing common-bin density estimates."""
    config, blocks = _fixture()
    original, _, _ = fit_fixed_data(config, blocks, 5.0, 15.0)
    edges = np.asarray(config["s_bin_edges"])
    extended = config | {"s_bin_edges": np.r_[edges, edges[-1] + np.diff(edges)[0]]}
    revised, _, _ = fit_fixed_data(extended, blocks, 5.0, 15.0)
    np.testing.assert_allclose(
        np.array(original["pmf_kcal_mol"], float),
        np.array(revised["pmf_kcal_mol"][:-1], float),
        rtol=0,
        atol=1e-12,
        equal_nan=True,
    )
    assert revised["frames_per_bin"][-1] == 0
