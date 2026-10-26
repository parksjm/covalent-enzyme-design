"""Regression checks for supervised path startup from a prior stopped state."""

import pytest

from scripts import run_nyl12_path_inspection as inspection_runner
from scripts.run_nyl12_path_inspection import _terminal_ready


def test_supervised_resume_waits_for_child_state_transition() -> None:
    """A stale stopped state must not terminate a newly launched resume child."""
    assert not _terminal_ready("stopped", True, False, False)
    assert not _terminal_ready("running", True, True, False)
    assert _terminal_ready("stopped", True, True, False)
    assert _terminal_ready("completed", True, True, True)
    assert _terminal_ready("stopped", True, False, True)
    assert _terminal_ready("stopped", False, False, False)


def test_secondary_publisher_delay_is_explicit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    delays: list[float] = []
    monkeypatch.setattr(inspection_runner.time, "sleep", delays.append)
    inspection_runner._pause_before_publish(25)
    assert delays == [25]
    with pytest.raises(ValueError, match="negative"):
        inspection_runner._pause_before_publish(-1)
