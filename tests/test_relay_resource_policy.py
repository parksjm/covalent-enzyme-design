"""Resource selection switches only for gain exceeding declared drift bounds."""

from __future__ import annotations
import pytest
from scripts.finalize_relay_affinity import choose_affinity


@pytest.mark.parametrize(
    "before,physical,after,winner",
    [
        (100, 106, 101, "physical_28"),
        (100, 102, 101, "all_56_logical"),
        (100, 110, 108, "all_56_logical"),
        (100, 90, 100, "all_56_logical"),
    ],
)
def test_affinity_selection_respects_noise_and_drift(
    before: float, physical: float, after: float, winner: str
) -> None:
    assert choose_affinity(before, physical, after)["winner"] == winner


def test_affinity_selection_rejects_missing_throughput() -> None:
    with pytest.raises(ValueError, match="positive"):
        choose_affinity(0, 10, 10)
