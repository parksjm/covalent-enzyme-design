"""Tests for direct proton-transfer and attack coordinate driving."""

import pytest

from nyl12_covalent.coupled_ti1_drive import (
    attack_only_schedule,
    attack_with_proton_identity_restraints,
    asynchronous_extension_schedule,
    asynchronous_restraint,
    asynchronous_schedule,
    coupled_restraints,
    coupled_schedule,
    ti1_like_geometry,
)


def test_schedule_advances_both_coordinates_in_small_steps() -> None:
    rows = coupled_schedule(-1.09, 3.23)
    assert rows[-1]["proton_target_A"] == 0.75
    assert rows[-1]["attack_target_A"] == 1.55
    assert all(
        later["proton_target_A"] > earlier["proton_target_A"]
        and later["attack_target_A"] < earlier["attack_target_A"]
        for earlier, later in zip(rows, rows[1:])
    )
    starts = [{"proton_target_A": -1.09, "attack_target_A": 3.23}, *rows[:-1]]
    assert (
        max(
            later["proton_target_A"] - earlier["proton_target_A"]
            for earlier, later in zip(starts, rows)
        )
        <= 0.10
    )
    assert (
        max(
            earlier["attack_target_A"] - later["attack_target_A"]
            for earlier, later in zip(starts, rows)
        )
        <= 0.10
    )


def test_restraints_include_only_two_reaction_coordinates() -> None:
    text = coupled_restraints(0.25, 1.8)
    assert text.count("&rst") == 2
    assert "iat=13916,13917,13905,13917" in text
    assert "iat=13916,19916,0,0,0,0,0,0" in text
    assert "restraintmask" not in text


def test_asynchronous_schedule_advances_one_composite_coordinate() -> None:
    rows = asynchronous_schedule(-1.09, 3.23)
    assert rows[-1]["composite_target_A"] == pytest.approx(-0.80)
    starts = [-1.09 - 3.23, *[row["composite_target_A"] for row in rows[:-1]]]
    assert (
        max(row["composite_target_A"] - start for start, row in zip(starts, rows))
        <= 0.15
    )
    assert "proton_target_A" not in rows[0]
    assert "attack_target_A" not in rows[0]


def test_asynchronous_restraint_is_one_three_distance_collective_variable() -> None:
    text = asynchronous_restraint(-2.0)
    assert text.count("&rst") == 1
    assert "iat=13916,13917,13905,13917,13916,19916,0,0" in text
    assert "rstwt=1,-1,-1,0" in text


def test_asynchronous_extension_preserves_spacing_and_labels() -> None:
    rows = asynchronous_extension_schedule(-0.8, 0.15, 25, maximum_windows=3)
    assert [row["label"] for row in rows] == [
        "asynchronous_25",
        "asynchronous_26",
        "asynchronous_27",
    ]
    assert [row["composite_target_A"] for row in rows] == pytest.approx(
        [-0.65, -0.50, -0.35]
    )


@pytest.mark.parametrize(
    ("increment", "first", "windows"),
    [(0.0, 25, 3), (0.15, 0, 3), (0.15, 25, 0)],
)
def test_asynchronous_extension_rejects_invalid_limits(
    increment: float, first: int, windows: int
) -> None:
    with pytest.raises(ValueError):
        asynchronous_extension_schedule(-0.8, increment, first, windows)


def test_ti1_like_geometry_requires_bond_and_oxyanion_rehybridization() -> None:
    geometry = {
        "OG1_carbonyl_C_A": 1.72,
        "N_HG1_A": 1.03,
        "OG1_HG1_A": 2.20,
        "carbonyl_C_O_A": 1.33,
    }
    assert ti1_like_geometry(geometry)
    geometry["OG1_carbonyl_C_A"] = 2.10
    assert not ti1_like_geometry(geometry)


def test_attack_only_schedule_uses_small_sequential_steps() -> None:
    rows = attack_only_schedule(2.1036)
    assert rows[-1]["attack_target_A"] == pytest.approx(1.50)
    starts = [2.1036, *[row["attack_target_A"] for row in rows[:-1]]]
    assert (
        max(start - row["attack_target_A"] for start, row in zip(starts, rows)) <= 0.075
    )


def test_attack_restraints_gate_all_terminal_nitrogen_protons() -> None:
    text = attack_with_proton_identity_restraints(1.75)
    assert text.count("&rst") == 4
    assert "iat=13916,19916" in text
    assert "iat=13905,13917" in text
    assert "iat=13905,13906" in text
    assert "iat=13905,13907" in text
    assert text.count("r2=0.9000000000, r3=1.2000000000") == 3
