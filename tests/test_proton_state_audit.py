"""Tests for cpptraj proton-state distance summaries."""

import json
from pathlib import Path

import pytest
from pytest import approx

from nyl12_covalent.proton_state_audit import (
    REQUIRED_COLUMNS,
    classify_hydrogen,
    summarize_stage,
    summarize_tables,
    write_summary,
)


def _write_table(path: Path, rows: list[str]) -> Path:
    path.write_text(
        " ".join(REQUIRED_COLUMNS) + "\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def test_classifies_boundary_values_and_ambiguous() -> None:
    assert classify_hydrogen(1.25, 1.30) == "N-bound"
    assert classify_hydrogen(1.30, 1.25) == "OG1-bound"
    assert classify_hydrogen(1.27, 1.27) == "ambiguous"


def test_summarizes_states_intended_frames_and_distances(tmp_path: Path) -> None:
    path = _write_table(
        tmp_path / "stage_a.dat",
        [
            "1 1.0 1.5 1.1 1.6 1.0 1.7 2.4 1.2",
            "2 1.4 1.1 1.0 1.5 1.0 1.5 2.6 1.4",
        ],
    )
    summary = summarize_stage(path)
    assert summary["frame_count"] == 2
    assert summary["state_counts"]["HG1"] == {"N-bound": 1, "OG1-bound": 1}
    assert summary["intended_state_frames"] == 1
    assert summary["OG1_C_le_2_5_frames"] == 1
    stats = summary["distance_stats"]["C_O"]
    assert stats["min"] == 1.2
    assert stats["mean"] == approx(1.3)
    assert stats["max"] == 1.4


def test_writes_multiple_stage_json(tmp_path: Path) -> None:
    first = _write_table(tmp_path / "pre.dat", ["1 1 1.5 1 1.5 1 1.5 2.0 1"])
    second = _write_table(tmp_path / "post.dat", ["1 1 1.5 1 1.5 1 1.5 2.0 1"])
    output = tmp_path / "nested" / "summary.json"
    write_summary([first, second], output)
    assert [item["stage"] for item in json.loads(output.read_text())["stages"]] == [
        "pre",
        "post",
    ]


def test_rejects_wrong_header_and_empty_path_list(tmp_path: Path) -> None:
    bad = tmp_path / "bad.dat"
    bad.write_text("Frame other\n1 2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected header"):
        summarize_stage(bad)
    with pytest.raises(ValueError, match="at least one"):
        summarize_tables([])


def test_rejects_malformed_row(tmp_path: Path) -> None:
    path = _write_table(tmp_path / "bad.dat", ["1 2"])
    with pytest.raises(ValueError, match="expected 9 fields"):
        summarize_stage(path)
