"""Tests for stable fixed-grid QM/MM umbrella PMF analysis."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.relay_fixed_pmf import contiguous_bootstrap, fit_samples


TEMPERATURE = 300.0
EDGES = np.array([0.0, 1.0, 2.0, 3.0])
REFERENCE = 0.5


def identical_umbrella_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return two identical zero-bias umbrellas with known bin populations."""
    one_window = np.array([0.2] * 4 + [1.2] * 2 + [2.2] * 6)
    s = np.concatenate([one_window, one_window])
    counts = np.array([len(one_window), len(one_window)])
    reduced = np.zeros((2, len(s)))
    return reduced, counts, s


def test_identical_umbrellas_match_analytic_histogram_and_report_counts() -> None:
    reduced, counts, s = identical_umbrella_data()

    result = fit_samples(reduced, counts, s, EDGES, TEMPERATURE, REFERENCE)

    kbt = 0.00198720425864083 * TEMPERATURE
    expected = -kbt * np.log(np.array([8.0, 4.0, 12.0]) / 8.0)
    expected -= expected[0]
    np.testing.assert_allclose(result["pmf_kcal_mol"], expected, atol=1e-7)
    assert result["frames_per_bin"] == [8, 4, 12]
    assert result["N_k"] == counts.tolist()
    assert result["reference_s"] == REFERENCE
    assert result["pmf_kcal_mol"][0] == pytest.approx(0.0, abs=1e-12)
    assert result["sampling_converged"] is False
    assert result["global_barrier_accepted"] is False
    assert "no temporal-correlation correction" in result["weight_ess_scope"]


def test_missing_bin_is_null_and_has_zero_support_ess() -> None:
    reduced, counts, s = identical_umbrella_data()
    s = np.where(s == 2.2, 0.2, s)

    result = fit_samples(reduced, counts, s, EDGES, TEMPERATURE, REFERENCE)

    assert result["pmf_kcal_mol"][2] is None
    assert result["frames_per_bin"][2] == 0
    assert result["weight_ess_per_bin"][2] == 0


def test_common_reference_is_fixed_and_high_dynamic_common_offsets_are_stable() -> None:
    reduced, counts, s = identical_umbrella_data()
    baseline = fit_samples(reduced, counts, s, EDGES, TEMPERATURE, 1.5)
    # Frame-dependent reduced potentials span a range large enough to make
    # direct exponentiation underflow. Log-sum-exp should retain finite results.
    offsets = np.linspace(-10000.0, 10000.0, len(s))
    shifted = reduced + offsets[None, :]

    result = fit_samples(shifted, counts, s, EDGES, TEMPERATURE, 1.5)

    kbt = 0.00198720425864083 * TEMPERATURE
    bins = np.searchsorted(EDGES, s, side="right") - 1
    # These inputs are umbrella potentials, so the unrestrained target's
    # importance weight is proportional to exp(+u_bias).
    log_mass = np.array(
        [np.logaddexp.reduce(offsets[bins == i]) for i in range(len(EDGES) - 1)]
    )
    expected = -kbt * (log_mass - np.log(np.diff(EDGES)))
    expected -= expected[1]
    np.testing.assert_allclose(result["pmf_kcal_mol"], expected, atol=1e-7)
    assert result["pmf_kcal_mol"] != baseline["pmf_kcal_mol"]
    assert result["reference_s"] == baseline["reference_s"] == 1.5
    assert np.isfinite(result["global_weight_ess"])
    assert result["global_weight_ess"] > 0


@pytest.mark.parametrize(
    ("reduced", "counts", "s", "edges", "temperature", "reference_s"),
    [
        (np.zeros((2, 24)), np.array([10, 13]), np.zeros(24), EDGES, TEMPERATURE, 0.5),
        (
            np.zeros((2, 24)),
            np.array([12.0, 12.0]),
            np.zeros(24),
            EDGES,
            TEMPERATURE,
            0.5,
        ),
        (np.zeros((2, 23)), np.array([12, 12]), np.zeros(24), EDGES, TEMPERATURE, 0.5),
        (
            np.zeros((2, 24)),
            np.array([12, 12]),
            np.full(24, np.nan),
            EDGES,
            TEMPERATURE,
            0.5,
        ),
        (np.zeros((2, 24)), np.array([12, 12]), np.zeros(24), EDGES, 0.0, 0.5),
        (np.zeros((2, 24)), np.array([12, 12]), np.zeros(24), EDGES, TEMPERATURE, 0.0),
        (
            np.zeros((2, 24)),
            np.array([12, 12]),
            np.zeros(24),
            np.array([0, 2, 1]),
            TEMPERATURE,
            0.5,
        ),
        (
            np.full((2, 24), np.inf),
            np.array([12, 12]),
            np.zeros(24),
            EDGES,
            TEMPERATURE,
            0.5,
        ),
    ],
)
def test_fit_validation_errors(
    reduced: np.ndarray,
    counts: np.ndarray,
    s: np.ndarray,
    edges: np.ndarray,
    temperature: float,
    reference_s: float,
) -> None:
    with pytest.raises(ValueError):
        fit_samples(reduced, counts, s, edges, temperature, reference_s)


def test_edge_coverage_and_reference_center_are_validated() -> None:
    reduced, counts, s = identical_umbrella_data()
    with pytest.raises(ValueError, match="strictly cover"):
        fit_samples(
            reduced, counts, s, np.array([0.0, 1.0, 2.0, 2.2]), TEMPERATURE, REFERENCE
        )
    with pytest.raises(ValueError, match="exactly one"):
        fit_samples(reduced, counts, s, EDGES, TEMPERATURE, 1.0)


def test_bootstrap_is_deterministic_and_common_reference_stays_fixed() -> None:
    reduced, counts, s = identical_umbrella_data()
    kwargs = dict(
        reduced=reduced,
        counts=counts,
        s=s,
        edges=EDGES,
        temperature=TEMPERATURE,
        reference_s=0.5,
        chunk_frames=2,
        replicates=50,
        seed=179,
    )

    first = contiguous_bootstrap(**kwargs, workers=1)
    repeated = contiguous_bootstrap(**kwargs, workers=1)

    assert first == repeated
    assert first["support_fraction"][1] < 1.0
    assert first["lower_95_kcal_mol"][1] is None
    assert all(
        draw[0] is None or draw[0] == pytest.approx(0.0, abs=1e-12)
        for draw in first["replicates_kcal_mol"]
    )
    assert first["replicates"] == 50
    assert first["workers"] == 1
    assert first["interval_scope"].startswith("Pointwise conditional")


def test_bootstrap_two_workers_matches_one_worker() -> None:
    reduced, counts, s = identical_umbrella_data()
    kwargs = dict(
        reduced=reduced,
        counts=counts,
        s=s,
        edges=EDGES,
        temperature=TEMPERATURE,
        reference_s=0.5,
        chunk_frames=2,
        replicates=50,
        seed=251,
    )

    serial = contiguous_bootstrap(**kwargs, workers=1)
    parallel = contiguous_bootstrap(**kwargs, workers=2)

    assert serial["replicates_kcal_mol"] == parallel["replicates_kcal_mol"]
    assert serial["support_fraction"] == parallel["support_fraction"]
    assert parallel["workers"] == 2


def test_bootstrap_full_support_produces_intervals() -> None:
    one_window = np.array([0.2, 1.2, 2.2] * 4)
    s = np.concatenate([one_window, one_window])
    counts = np.array([12, 12])

    result = contiguous_bootstrap(
        np.zeros((2, len(s))),
        counts,
        s,
        EDGES,
        TEMPERATURE,
        REFERENCE,
        chunk_frames=3,
        replicates=50,
        workers=1,
        seed=47,
    )

    assert result["support_fraction"] == [1.0, 1.0, 1.0]
    assert all(value is not None for value in result["lower_95_kcal_mol"])
    assert all(value is not None for value in result["upper_95_kcal_mol"])


def test_disconnected_windows_are_null_in_point_and_bootstrap() -> None:
    one_window = np.array([0.2] * 12)
    other_window = np.array([2.2] * 12)
    s = np.concatenate([one_window, other_window])
    counts = np.array([12, 12])
    reduced = np.array(
        [
            np.concatenate([np.zeros(12), np.full(12, 1000.0)]),
            np.concatenate([np.full(12, 1000.0), np.zeros(12)]),
        ]
    )

    point = fit_samples(reduced, counts, s, EDGES, TEMPERATURE, REFERENCE)
    boot = contiguous_bootstrap(
        reduced,
        counts,
        s,
        EDGES,
        TEMPERATURE,
        REFERENCE,
        chunk_frames=2,
        replicates=50,
        workers=1,
        seed=59,
    )

    assert point["connected"] is False
    assert point["pmf_kcal_mol"] == [None, None, None]
    assert boot["disconnected_replicates"] == 50
    assert boot["support_fraction"] == [0.0, 0.0, 0.0]
    assert boot["lower_95_kcal_mol"] == [None, None, None]


@pytest.mark.parametrize(
    ("chunk_frames", "replicates", "workers"),
    [(0, 50, 1), (13, 50, 1), (7, 50, 1), (2, 49, 1), (2, 50, 0), (2, 50, 9)],
)
def test_bootstrap_validation_errors(
    chunk_frames: int, replicates: int, workers: int
) -> None:
    reduced, counts, s = identical_umbrella_data()

    with pytest.raises(ValueError):
        contiguous_bootstrap(
            reduced,
            counts,
            s,
            EDGES,
            TEMPERATURE,
            REFERENCE,
            chunk_frames=chunk_frames,
            replicates=replicates,
            workers=workers,
            seed=1,
        )
