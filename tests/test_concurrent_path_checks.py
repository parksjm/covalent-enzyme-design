"""Tests for mapped gates against composite-coordinate degeneracy."""

import pytest
from nyl12_covalent.concurrent_path_checks import candidate_geometry, chemistry_retained


def geometry() -> dict[str, float]:
    """Return a plausible mapped candidate."""
    return {
        "N_H1_A": 1.04,
        "N_H2_A": 1.04,
        "OG1_H1_A": 2.1,
        "OG1_H2_A": 2.1,
        "N_HG1_A": 1.05,
        "OG1_HG1_A": 1.8,
        "OG1_carbonyl_C_A": 1.65,
        "carbonyl_C_O_A": 1.33,
        "carbonyl_C_N_A": 1.45,
    }


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("N_H1_A", 1.4),
        ("N_H2_A", 1.4),
        ("OG1_H1_A", 1.0),
        ("OG1_H2_A", 1.0),
        ("OG1_carbonyl_C_A", 1.1),
        ("carbonyl_C_O_A", 1.7),
        ("carbonyl_C_N_A", 1.9),
        ("N_HG1_A", 0.6),
    ],
)
def test_reject_components(key: str, value: float) -> None:
    values = geometry()
    values[key] = value
    assert not chemistry_retained(values)
    assert not candidate_geometry(values)


def test_shared_proton_allowed_but_not_candidate() -> None:
    values = geometry()
    values.update(N_HG1_A=1.35, OG1_HG1_A=1.35)
    assert chemistry_retained(values)
    assert not candidate_geometry(values)
    values.update(N_HG1_A=1.8, OG1_HG1_A=1.8)
    assert not chemistry_retained(values)


def test_candidate_requires_attack_and_oxyanion() -> None:
    values = geometry()
    assert candidate_geometry(values)
    values["OG1_carbonyl_C_A"] = 2.1
    assert not candidate_geometry(values)
    values["OG1_carbonyl_C_A"] = 1.65
    values["carbonyl_C_O_A"] = 1.23
    assert not candidate_geometry(values)
