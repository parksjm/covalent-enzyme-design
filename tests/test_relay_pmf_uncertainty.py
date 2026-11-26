"""Verify pointwise errors against independent multinomial statistics."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent import relay_pmf_uncertainty as module
from nyl12_covalent.relay_pmf_uncertainty import histogram_uncertainty


def independent_histogram(repeats: int = 1) -> dict:
    """Use identical sampled states with exactly known histogram covariance."""
    s = np.repeat([0.5, 1.5, 3.5], np.array([30, 20, 50]) * repeats)
    return histogram_uncertainty(
        np.zeros((2, len(s))),
        np.array([50, 50]) * repeats,
        s,
        np.arange(5),
        300.0,
        0.5,
    )


def test_covariance_agrees_with_multinomial_and_preserves_gap() -> None:
    result = independent_histogram()
    kbt = module.KB * 300
    assert result["standard_error_kcal_mol"][0] == pytest.approx(0)
    assert result["standard_error_kcal_mol"][1] == pytest.approx(
        kbt * np.sqrt(1 / 20 + 1 / 30), rel=1e-7
    )
    assert result["pmf_kcal_mol"][1] == pytest.approx(-kbt * np.log(20 / 30))
    for key in ("pmf_kcal_mol", "standard_error_kcal_mol", "lower_95_kcal_mol"):
        assert result[key][2] is None
    assert result["selected_frames_per_bin"] == [30, 20, 0, 50]
    assert result["weight_ess_per_bin"] == pytest.approx([30, 20, 0, 50])
    assert result["contributing_windows_per_bin"] == [1, 1, 0, 1]
    assert not result["sampling_converged"]


def test_standard_error_scales_as_inverse_square_root_of_sample_size() -> None:
    first, second = independent_histogram(), independent_histogram(4)
    np.testing.assert_allclose(
        np.array(first["standard_error_kcal_mol"], float),
        np.array(second["standard_error_kcal_mol"], float) * 2,
        equal_nan=True,
    )


def test_nonuniform_bins_report_density_not_probability_mass() -> None:
    result = histogram_uncertainty(
        np.zeros((2, 100)),
        np.array([50, 50]),
        np.repeat([0.5, 2], 50),
        np.array([0, 1, 3]),
        300,
        0.5,
    )
    assert result["pmf_kcal_mol"][1] == pytest.approx(module.KB * 300 * np.log(2))


@pytest.mark.parametrize(
    ("counts", "edges", "temperature", "reference", "message"),
    [
        ([0, 4], [0, 1, 2], 300, 0.5, "positive counts"),
        ([2, 2], [0, 1, 2], -300, 0.5, "temperature"),
        ([2, 2], [0, 2, 1], 300, 0.5, "edges"),
        ([2, 2], [0, 1, 2], 300, 9, "Reference"),
    ],
)
def test_invalid_inputs_fail(
    counts: list[int],
    edges: list[int],
    temperature: float,
    reference: float,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        histogram_uncertainty(
            np.zeros((2, 4)),
            np.array(counts),
            np.repeat(0.5, 4),
            np.array(edges),
            temperature,
            reference,
        )


def test_same_source_selection_is_enforced(monkeypatch: pytest.MonkeyPatch) -> None:
    summary = dict(
        selected_frames=2,
        raw_frames=4,
        discarded_frames=0,
        statistical_inefficiency=2.0,
    )
    monkeypatch.setattr(
        module, "select_uncorrelated", lambda _: (np.array([0, 2]), summary)
    )
    blocks = [dict(window="w0", values=np.zeros((4, 3)))]
    report = dict(N_k=[2], blocks=[dict(window="w0", **summary)])
    samples, counts = module.selected_samples([dict(label="w0")], blocks, report)
    assert samples.shape == (2, 3)
    assert counts.tolist() == [2]
    with pytest.raises(ValueError, match="summaries"):
        module.selected_samples([dict(label="w0")], blocks, dict(blocks=[]))
    report["N_k"] = [3]
    with pytest.raises(ValueError, match="counts"):
        module.selected_samples([dict(label="w0")], blocks, report)
    report["blocks"][0]["discarded_frames"] = 1
    with pytest.raises(ValueError, match="selection"):
        module.selected_samples([dict(label="w0")], blocks, report)


def test_report_identity_and_projection_guards(monkeypatch: pytest.MonkeyPatch) -> None:
    bias = dict(
        feature_names=["q_acyl", "q_act"],
        centers_A=[0, 0],
        force_constants=[1, 1],
        outer_extent_A=10,
    )
    samples = np.zeros((4, 3))
    samples[:, 0] = [0.5, 0.5, 1.5, 1.5]
    config = dict(
        temperature_K=300,
        windows=[dict(bias=bias)],
        path_3cv_A=[[0, 0, 0], [2, 0, 0]],
        s_bin_edges=[0, 0.5, 1],
    )
    report = dict(
        temperature_K=300,
        pmf=dict(pmf_kcal_mol=[0, 0]),
        projection_diagnostics=dict(s=[0.25, 0.25, 0.75, 0.75]),
        checkpoint=dict(plot_common_es_reference_s=0.25),
    )
    monkeypatch.setattr(
        module, "selected_samples", lambda *args: (samples, np.array([4]))
    )
    monkeypatch.setattr(
        module,
        "histogram_uncertainty",
        lambda *args: dict(centers=[0.25, 0.75], pmf_kcal_mol=[0, 0]),
    )
    assert module.report_uncertainty(config, [], report)[
        "source_estimator_identity_verified"
    ]
    report["pmf"]["pmf_kcal_mol"] = [0, 2]
    with pytest.raises(ValueError, match="identity"):
        module.report_uncertainty(config, [], report)
    report["projection_diagnostics"]["s"][0] = -5
    with pytest.raises(ValueError, match="projections"):
        module.report_uncertainty(config, [], report)
    config["temperature_K"] = 301
    with pytest.raises(ValueError, match="temperature"):
        module.report_uncertainty(config, [], report)
