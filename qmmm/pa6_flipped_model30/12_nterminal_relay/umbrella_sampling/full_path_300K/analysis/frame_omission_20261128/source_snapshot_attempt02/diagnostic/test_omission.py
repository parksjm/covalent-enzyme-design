"""Small independent scientific qualification for the deletion diagnostic."""

import numpy as np
import pytest
from scipy.special import logsumexp

from omission_core import exact_refit, omit_observation, profile


def test_exactly_one_provenance_match_and_state_decrement() -> None:
    q = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [1.0, 2.0, 3.0], [7.0, 8.0, 9.0]])
    keep, counts, index = omit_observation(q, np.array([2, 2]), 1, q[0])
    assert index == 2 and keep.tolist() == [True, True, False, True]
    assert counts.tolist() == [2, 1]
    with pytest.raises(ValueError, match="exactly one"):
        omit_observation(q, np.array([2, 2]), 1, np.zeros(3))
    with pytest.raises(ValueError, match="counts"):
        omit_observation(q, np.array([1, 2]), 1, q[0])
    with pytest.raises(ValueError, match="sampled"):
        omit_observation(q, np.array([3, 1]), 1, q[-1])


def test_fixed_weight_bin_shift_matches_analytical_removal() -> None:
    weights = np.array([0.2, 0.05, 0.6, 0.15])
    s, edges = np.array([0.1, 0.1, 0.3, 0.3]), np.array([0.0, 0.2, 0.4])
    base = profile(s, np.log(weights), edges, 0, 0.596, np.arange(4))
    keep = np.array([True, True, False, True])
    omitted = profile(
        s[keep], np.log(weights[keep]), edges, 0, 0.596, np.arange(4)[keep]
    )
    np.testing.assert_allclose(omitted["pmf"][1] - base["pmf"][1], -0.596 * np.log(0.2))
    assert omitted["pmf"][0] == 0
    np.testing.assert_allclose(omitted["mass"].sum(), 1)
    assert omitted["dominant_original_index"][1] == 3
    with pytest.raises(ValueError, match="reference"):
        profile(s[2:], np.log(weights[2:]), edges, 0, 0.596, np.arange(2))
    with pytest.raises(ValueError, match="inputs"):
        profile(s, np.log(weights), edges, 2, 0.596, np.arange(4))


def test_synthetic_unbiased_states_recover_uniform_observation_weights() -> None:
    # Both states describe exactly the same unbiased Hamiltonian with known
    # constant energy offset. Their unequal sample counts must be honored.
    counts = np.array([3, 2])
    u = np.tile(np.array([0.0, 1.7])[:, None], (1, 5))
    result = exact_refit(u, counts, np.array([0.0, 1.7]))
    np.testing.assert_allclose(result["f_k"], [0.0, 1.7], atol=1e-10)
    np.testing.assert_allclose(result["logweights"], -np.log(5), atol=1e-10)
    assert result["connected_components_at_1e_minus12"] == 1
    assert result["fixed_point_residual_kBT"] < 1e-10
    np.testing.assert_allclose(logsumexp(result["logweights"]), 0, atol=1e-14)
    with pytest.raises(ValueError, match="inputs"):
        exact_refit(u, np.array([3, 1]), np.zeros(2))
