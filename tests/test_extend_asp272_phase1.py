"""Focused safeguards for continuation of the original phase I path."""

import json
from pathlib import Path

import pytest

from provenance.asp272_relay_20261025.extend_phase1 import schedule, validate_status


def test_schedule_appends_two_original_interval_steps() -> None:
    windows = schedule(3.1410380093, -0.9988674481)

    assert [item["label"] for item in windows] == ["stage1_17", "stage1_18"]
    assert windows[0]["attack_target_A"] == pytest.approx(1.4505601244)
    assert windows[1]["attack_target_A"] == pytest.approx(1.3511202488)
    assert windows[0]["q1_target_A"] == pytest.approx(0.912429215)
    assert windows[1]["q1_target_A"] == pytest.approx(1.02485843)
    assert windows[1]["attack_target_A"] - windows[0][
        "attack_target_A"
    ] == pytest.approx((1.55 - 3.1410380093) / 16)
    assert windows[0]["attack_target_A"] - 1.55 == pytest.approx(
        (1.55 - 3.1410380093) / 16
    )
    assert windows[1]["q1_target_A"] - windows[0]["q1_target_A"] == pytest.approx(
        (0.8 + 0.9988674481) / 16
    )


def test_status_validation_rejects_incomplete_phase(tmp_path: Path) -> None:
    (tmp_path / "status.json").write_text(json.dumps({"phase1_complete": False}))

    with pytest.raises(ValueError, match="Phase I must be complete and stopped"):
        validate_status(tmp_path)
