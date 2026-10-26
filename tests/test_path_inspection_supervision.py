"""Regression checks for supervised path startup from a prior stopped state."""

from scripts.run_nyl12_path_inspection import _terminal_ready


def test_supervised_resume_waits_for_child_state_transition() -> None:
    """A stale stopped state must not terminate a newly launched resume child."""
    assert not _terminal_ready("stopped", True, False, False)
    assert not _terminal_ready("running", True, True, False)
    assert _terminal_ready("stopped", True, True, False)
    assert _terminal_ready("completed", True, True, True)
    assert _terminal_ready("stopped", True, False, True)
    assert _terminal_ready("stopped", False, False, False)
