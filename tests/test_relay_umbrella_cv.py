"""Focused checks for relay CV measurement and generalized Amber biases."""

from __future__ import annotations

import numpy as np
import pytest

from nyl12_covalent.neb_umbrella_bias import HarmonicWindow, bias_energy
from nyl12_covalent.relay_umbrella_cv import (
    DEFAULT_BIAS_FEATURES,
    FEATURES,
    DistanceDifference,
    RelayBiasWindow,
    amber_bias_energy,
    cross_bias_energies,
    iter_raw_features,
    raw_features,
)


def _relay_frame() -> np.ndarray:
    xyz = np.zeros((19929, 3), dtype=float)
    xyz[19916 - 1] = (0.0, 0.0, 0.0)
    xyz[19929 - 1] = (2.0, 0.0, 0.0)
    xyz[13916 - 1] = (1.0, 0.0, 0.0)
    xyz[13905 - 1] = (0.0, 3.0, 0.0)
    xyz[13917 - 1] = (0.0, 2.0, 0.0)
    return xyz


def test_distance_difference_matches_analytic_geometry() -> None:
    xyz = np.array([[0, 0, 0], [3, 0, 0], [0, 4, 0], [0, 0, 1]], dtype=float)
    feature = DistanceDifference("analytic", (1, 2, 3, 4))
    assert feature.value(xyz) == pytest.approx(3.0 - np.sqrt(17.0))


def test_raw_features_include_all_three_relay_coordinates_and_q_h() -> None:
    values = raw_features(_relay_frame())
    assert values["q_acyl"] == pytest.approx(1.0)
    assert values["q_act"] == pytest.approx(np.sqrt(5.0) - 1.0)
    assert values["q_del"] == pytest.approx(1.0 - np.sqrt(8.0))
    assert values["q_h"] == pytest.approx(np.sqrt(5.0) - np.sqrt(8.0))
    assert tuple(name for name in values if name != "q_h") == (
        "q_acyl",
        "q_act",
        "q_del",
    )


def test_frame_iterator_emits_all_raw_features_for_each_frame() -> None:
    frames = iter_raw_features((_relay_frame(), _relay_frame()))
    records = list(frames)
    assert len(records) == 2
    assert all(set(record) == {"q_acyl", "q_act", "q_del", "q_h"} for record in records)


def test_triclinic_minimum_image_and_coordinate_rewrapping() -> None:
    box = (10.0, 12.0, 14.0, 90.0, 90.0, 60.0)
    original = _relay_frame()
    wrapped_equivalent = original.copy()
    wrapped_equivalent[13917 - 1] += np.array([10.0, 0.0, 0.0])
    first = raw_features(original, box)
    second = raw_features(wrapped_equivalent, box)
    for name in ("q_acyl", "q_act", "q_del", "q_h"):
        assert second[name] == pytest.approx(first[name])
    assert first["q_h"] == pytest.approx(first["q_act"] + first["q_del"])


def test_periodic_cell_validation_and_frame_box_counts() -> None:
    with pytest.raises(ValueError, match="six finite"):
        raw_features(_relay_frame(), (10.0, 10.0))
    with pytest.raises(ValueError):
        list(iter_raw_features((_relay_frame(),), ()))


def test_native_restraint_atom_roles_and_default_q_h_compatibility() -> None:
    window = RelayBiasWindow(
        ("q_acyl", "q_act"), (-0.25, 1.10), (25.0, 25.0)
    )
    text = window.restraint_text()
    assert text.count("&rst") == 2
    assert "iat=19916,19929,19916,13916,0,0,0,0" in text
    assert "iat=13916,13917,13905,13917,0,0,0,0" in text
    assert text.count("rstwt=1,-1,0,0") == 2
    assert "rk2=25, rk3=25" in text

    legacy = RelayBiasWindow()
    assert legacy.feature_names == DEFAULT_BIAS_FEATURES == ("q_acyl", "q_h")
    assert "iat=13916,13917,19929,13917,0,0,0,0" in legacy.restraint_text()


def test_native_amber_energy_has_no_half_factor_and_keeps_tails() -> None:
    assert amber_bias_energy(0.2, 0.0, 25.0, 10.0) == pytest.approx(1.0)
    assert amber_bias_energy(11.0, 0.0, 25.0, 10.0) == pytest.approx(3000.0)


def test_cross_bias_maps_arbitrary_chart_axes_for_every_frame() -> None:
    windows = (
        RelayBiasWindow(("q_acyl", "q_act"), (0.0, 1.0), (25.0, 25.0)),
        RelayBiasWindow(("q_acyl", "q_del"), (0.5, -1.0), (10.0, 20.0)),
    )
    features = {
        "q_acyl": np.array([0.0, 1.0]),
        "q_act": np.array([1.0, 2.0]),
        "q_del": np.array([-1.0, -2.0]),
    }
    matrix = cross_bias_energies(features, windows)
    assert matrix.shape == (2, 2)
    assert matrix[0].tolist() == pytest.approx([0.0, 50.0])
    assert matrix[1].tolist() == pytest.approx([2.5, 22.5])


def test_legacy_default_bias_matches_existing_two_cv_helper() -> None:
    old = HarmonicWindow(-0.5, 0.25, 25.0, 25.0)
    new = RelayBiasWindow(centers_A=(-0.5, 0.25), force_constants=(25.0, 25.0))
    features = {"q_acyl": np.array([-0.2]), "q_h": np.array([0.55])}
    observed = cross_bias_energies(features, (new,))[0, 0]
    assert observed == pytest.approx(bias_energy((-0.2, 0.55), old))


def test_invalid_feature_ids_and_geometry_fail() -> None:
    with pytest.raises(ValueError):
        DistanceDifference("bad", (1, 0, 2, 3))
    with pytest.raises(ValueError):
        DistanceDifference("bad", (1, 2, 3, 5)).value(np.zeros((4, 3)))
    with pytest.raises(ValueError):
        DistanceDifference("bad", (1, 2, 3, 4)).value(np.zeros((4,)))
    coordinates = np.zeros((4, 3))
    coordinates[3, 0] = np.nan
    with pytest.raises(ValueError):
        DistanceDifference("bad", (1, 2, 3, 4)).value(coordinates)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"feature_names": ("q_act", "q_act")},
        {"feature_names": ("q_act", "missing")},
        {"force_constants": (25.0, 0.0)},
        {"centers_A": (float("nan"), 0.0)},
    ],
)
def test_invalid_window_fails(kwargs: dict) -> None:
    values = {
        "feature_names": ("q_acyl", "q_act"),
        "centers_A": (0.0, 0.0),
        "force_constants": (25.0, 25.0),
    }
    values.update(kwargs)
    with pytest.raises(ValueError):
        RelayBiasWindow(**values)


def test_cross_bias_rejects_missing_or_invalid_feature_series() -> None:
    window = RelayBiasWindow(("q_acyl", "q_del"), (0.0, 0.0), (25.0, 25.0))
    with pytest.raises(ValueError, match="Missing frame features"):
        cross_bias_energies({"q_acyl": [0.0]}, (window,))
    with pytest.raises(ValueError, match="same-length"):
        cross_bias_energies({"q_acyl": [0.0, 1.0], "q_del": [0.0]}, (window,))
    with pytest.raises(ValueError, match="finite"):
        cross_bias_energies(
            {"q_acyl": [0.0], "q_del": [np.inf]}, (window,)
        )
    with pytest.raises(ValueError, match="At least one"):
        cross_bias_energies({}, ())


def test_registry_defines_exact_three_raw_features() -> None:
    assert tuple(name for name in FEATURES if name != "q_h") == (
        "q_acyl",
        "q_act",
        "q_del",
    )
