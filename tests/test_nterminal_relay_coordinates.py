"""Tests for relay coordinate, restraint and endpoint contracts."""

import numpy as np
import pytest

from nyl12_covalent import nterminal_relay_coordinates as relay
from nyl12_covalent.relay_release import ROLES


def coordinates() -> np.ndarray:
    """Build a full-sized coordinate array with distinct catalytic centers."""
    xyz = np.zeros((104922, 3))
    positions = {
        "C": (0, 0, 0),
        "O": (0, 1.2, 0),
        "OG": (1.5, 0, 0),
        "NL": (0, 0, 2.5),
        "CA": (-1.5, 0, 0),
        "N": (3, 0, 0),
        "H1": (3, 1, 0),
        "H2": (3, -1, 0),
        "HG": (0, 0, 2),
        "HN": (0, 0, 3.5),
    }
    for key, value in positions.items():
        xyz[ROLES[key] - 1] = value
    return xyz


def test_coordinate_matches_manual_expression() -> None:
    xyz = coordinates()
    for branch, acceptor in (("reverse_activation", "OG"), ("forward_acylation", "NL")):
        h = xyz[ROLES["HG"] - 1]
        a = xyz[ROLES[acceptor] - 1]
        expected = (
            np.linalg.norm(xyz[ROLES["N"] - 1] - h)
            - np.linalg.norm(a - h)
            + np.linalg.norm(xyz[ROLES["C"] - 1] - a)
        )
        assert relay.coordinate(xyz, branch) == pytest.approx(expected)


def test_restraint_and_explicit_hydrogen() -> None:
    text = relay.amber_restraint(2.0, "forward_acylation", 13906, 50)
    assert "iat=13905,13906,19929,13906,19916,19929,0,0" in text
    assert "rstwt=1,-1,1,0" in text
    assert "r1=0.0000000000" in text and "r4=4.0000000000" in text
    assert "rk2=50.0000000000" in text
    assert np.isfinite(relay.coordinate(coordinates(), "reverse_activation", 13907))


def test_schedule_bounds() -> None:
    rows = relay.schedule(1.0, 1.31, 0.15)
    assert [r["window"] for r in rows] == ["window_001", "window_002", "window_003"]
    values = [1.0] + [r["target"] for r in rows]
    assert values[-1] == pytest.approx(1.31)
    assert max(np.diff(values)) <= 0.15


@pytest.mark.parametrize(
    "args", [(1, 1, 0.1), (2, 1, 0.1), (1, 2, 0), (float("nan"), 2, 0.1)]
)
def test_invalid_schedule(args: tuple[float, float, float]) -> None:
    with pytest.raises(ValueError):
        relay.schedule(*args)


def test_input_validation() -> None:
    with pytest.raises(ValueError):
        relay.coordinate(np.zeros((3, 3)), "reverse_activation")
    with pytest.raises(ValueError):
        relay.coordinate(coordinates(), "other")
    with pytest.raises(ValueError):
        relay.coordinate(coordinates(), "reverse_activation", 19937)
    with pytest.raises(ValueError):
        relay.amber_restraint(1, "forward_acylation", force=0)
    xyz = coordinates()
    xyz[0, 0] = np.nan
    with pytest.raises(ValueError):
        relay.coordinate(xyz, "forward_acylation")
    with pytest.raises(ValueError):
        relay.candidate({}, "other")


def test_forward_ownership_and_identity_independence() -> None:
    g = relay.geometry(coordinates())
    assert relay.candidate(g, "forward_acylation")
    assert not relay.candidate(g, "reverse_activation")
    assert g["distances_A"]["N_NL"] > 0
    g["ownership"]["H1"], g["ownership"]["HG"] = (
        g["ownership"]["HG"],
        g["ownership"]["H1"],
    )
    assert relay.candidate(g, "forward_acylation")
    g["ownership"]["HG"] = ["N", "OG"]
    assert not relay.candidate(g, "forward_acylation")


def test_reverse_gates() -> None:
    g = relay.geometry(coordinates())
    g["distances_A"].update(C_OG=2.5, C_NL=1.4, C_O=1.2)
    g["ownership"].update(H1=["OG"], H2=["N"], HG=["N"], HN=["NL"])
    assert relay.candidate(g, "reverse_activation")
    g["ownership"]["HN"] = []
    assert not relay.candidate(g, "reverse_activation")
    g["ownership"]["HN"] = ["NL"]
    g["distances_A"]["C_O"] = 1.4
    assert not relay.candidate(g, "reverse_activation")
