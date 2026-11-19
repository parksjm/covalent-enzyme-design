"""Tests for the queued GPU launcher."""

import pytest

from scripts.run_when_gpu_idle import parse_pids


def test_parse_pids() -> None:
    assert parse_pids("123\n456\n") == [123, 456]
    assert parse_pids("N/A\n\n") == []
    with pytest.raises(ValueError):
        parse_pids("not-a-pid\n")
