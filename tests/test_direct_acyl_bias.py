"""Check Amber's quadratic and linear flat-bottom guard regions."""

import pytest

from nyl12_covalent.direct_acyl_bias import proton_guard_bias
from nyl12_covalent.direct_acyl_drive import direct_restraint


@pytest.mark.parametrize(
    "distance,expected",
    [(0.7, 2.625), (0.9, 0.125), (1.1, 0), (1.4, 1.125), (1.6, 5.625)],
)
def test_guard_regions(distance: float, expected: float) -> None:
    assert proton_guard_bias(
        direct_restraint(0.3, substrate_proton_guard=True), distance
    ) == pytest.approx(expected)


def test_guard_requires_audited_protocol() -> None:
    assert proton_guard_bias(direct_restraint(0.3), 1.4) is None
    with pytest.raises(ValueError, match="protocol"):
        proton_guard_bias(
            direct_restraint(0.3, substrate_proton_guard=True).replace(
                "rk2=50", "rk2=5"
            ),
            1.4,
        )
    with pytest.raises(ValueError, match="exactly"):
        proton_guard_bias(
            direct_restraint(0.3, substrate_proton_guard=True) + "&rst /", 1.4
        )
