"""Mathematical inventory and mass conversion checks for partial mode lifting."""

import numpy as np
import pytest

from nyl12_covalent.relay_ts_partial import lift_mode


def test_lift_preserves_full_dimension_and_mass_conversion() -> None:
    """A partial mode leaves every other full-shell coordinate at zero."""
    mode = lift_mode(
        np.array([1, 4]), np.array([0.6, 0.8]), np.array([1.0, 2.0, 1.0, 1.0, 0.5, 1.0])
    )
    expected = np.array([0.0, 1.2, 0.0, 0.0, 0.4, 0.0])
    np.testing.assert_allclose(mode, expected / np.linalg.norm(expected))


@pytest.mark.parametrize(
    "columns,vectors,masses",
    [
        ([1, 1], [1.0, 1.0], [1.0, 1.0]),
        ([2], [1.0], [1.0, 1.0]),
        ([-1], [1.0], [1.0, 1.0]),
        ([0], [0.0], [1.0, 1.0]),
        ([0], [np.nan], [1.0, 1.0]),
        ([0], [1.0], [-1.0, 1.0]),
        ([], [], [1.0, 1.0]),
    ],
)
def test_invalid_mode_fails(
    columns: list[int], vectors: list[float], masses: list[float]
) -> None:
    """Malformed or degenerate modes fail before any scientific evaluation."""
    with pytest.raises(ValueError):
        lift_mode(np.array(columns), np.array(vectors), np.array(masses))
