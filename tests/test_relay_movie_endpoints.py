"""Endpoint retention and bounded relay display-frame selection."""

import pytest

from scripts.export_neb_us_relay_release_vmd import select_optimizer_frames


@pytest.mark.parametrize("count", [1, 2, 46, 47, 48, 203, 348, 317])
def test_selection_retains_endpoints_within_movie_cap(count: int) -> None:
    selected = select_optimizer_frames(count)
    assert selected[0] == 0
    assert selected[-1] == count - 1
    assert len(selected) == min(count, 47)
    assert selected == sorted(set(selected))
    assert len(selected) + 1 <= 48


@pytest.mark.parametrize("count", [203, 348, 317])
def test_restart_reservation_preserves_native_endpoint(count: int) -> None:
    snapshot = select_optimizer_frames(count)
    selected = select_optimizer_frames(len(snapshot), 46)
    original = [snapshot[index] for index in selected]
    assert original[0] == 0
    assert original[-1] == count - 1
    assert len(original) + 2 == 48
    assert original == sorted(set(original))


@pytest.mark.parametrize(
    "count,limit", [(0, 47), (-1, 47), (2, 1), (1, 0), (1, -1), (True, 47), (2, 1.5)]
)
def test_selection_rejects_invalid_counts(count: int, limit: int) -> None:
    with pytest.raises(ValueError):
        select_optimizer_frames(count, limit)


def test_single_frame_needs_one_slot() -> None:
    assert select_optimizer_frames(1, 1) == [0]
