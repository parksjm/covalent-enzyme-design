"""Check matched WT unique-frame joins and rejection gates."""

import csv
import json
from pathlib import Path

import pytest

from nyl12_covalent.matched_wt_relays import STATES, join_row, summarize, validate_cases


def geometry() -> dict:
    return dict(d1=3.0, d2=3.0, d3=3.0, attack=100.0, projected_phi=-90.0)


def test_join_native_fail_and_identity() -> None:
    row = dict(
        case="x",
        frame=1,
        time_ps=510.0,
        water_count=1,
        internal_count=0,
        strict=True,
        **geometry(),
    )
    occupancy = {**row, "native_pass": "False", "tier": "A"}
    flags, category = join_row(row, occupancy)
    assert flags["strict_assisted_tier_ab"]
    assert not flags["strict_assisted_tier_ab_native"]
    assert category == "strict_joint_native_fail"
    for key, value in [
        ("time_ps", 511),
        ("case", "wrong"),
        ("frame", 2),
        ("d1", 3.1),
        ("native_pass", "bad"),
        ("tier", "B"),
    ]:
        with pytest.raises(ValueError):
            join_row(row, {**occupancy, key: value})


def test_complete_join_deduplicates_routes_and_repeated_screens(tmp_path: Path) -> None:
    occupancy = tmp_path / "occupancy"
    screen = tmp_path / "screen"
    occupancy.mkdir()
    screen.mkdir()
    cases, receipts = [], []
    for state in sorted(STATES):
        for replica in range(1, 5):
            name = f"{state}__replica_{replica:02d}"
            cases.append(
                dict(
                    name=name,
                    state=state,
                    topology="top",
                    topology_sha256="abc",
                    trajectories=["first", "second"],
                )
            )
            receipts.append(
                dict(
                    case=name,
                    total_frames=10000,
                    input_sha256={"top": "abc", "first": "def", "second": "ghi"},
                )
            )
            with (occupancy / f"{name}.csv").open("w", newline="") as handle:
                first = dict(
                    case=name,
                    frame=1,
                    time_ps=510.0,
                    tier="A",
                    native_pass="False",
                    **geometry(),
                )
                writer = csv.DictWriter(handle, fieldnames=list(first))
                writer.writeheader()
                for frame in range(1, 10001):
                    writer.writerow(
                        {**first, "frame": frame, "time_ps": 500 + frame * 10}
                    )
            with (screen / f"{name}_paths.jsonl").open("w") as handle:
                for frame in range(1, 10001, 5):
                    row = dict(
                        case=name,
                        frame=frame,
                        time_ps=500 + frame * 10,
                        strict=True,
                        water_count=1,
                        internal_count=0,
                        violation=0.0,
                        **geometry(),
                    )
                    handle.write((json.dumps(row) + "\n") * 2)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"cases": cases}))
    (screen / "summary.json").write_text(json.dumps({"cases": receipts}))
    for case in cases:
        (occupancy / f"{case['name']}_summary.json").write_text(
            json.dumps(
                {
                    "name": case["name"],
                    "state": case["state"],
                    "frames": 10000,
                    "input_sha256": {"top": "abc", "first": "def", "second": "ghi"},
                }
            )
        )
    summarize(config, occupancy, [screen, screen], tmp_path / "output")
    result = json.loads((tmp_path / "output/summary.json").read_text())
    assert len(result["replicas"]) == 12
    assert len(result["states"]) == 3
    assert all(row["sampled_frames"] == 2000 for row in result["replicas"])
    assert all(row["strict_assisted"] == 2000 for row in result["replicas"])
    assert all(row["strict_assisted_tier_ab_native"] == 0 for row in result["replicas"])
    assert all(
        row["candidate_class"] == "strict_joint_native_fail"
        for row in result["candidates"]
    )
    assert len(result["candidates"]) == 10
    stale_path = occupancy / f"{cases[0]['name']}_summary.json"
    stale = json.loads(stale_path.read_text())
    stale["input_sha256"]["second"] = "changed_water_coordinates_same_geometry"
    stale_path.write_text(json.dumps(stale))
    with pytest.raises(ValueError, match="MD input hash mismatch"):
        summarize(config, occupancy, [screen], tmp_path / "stale_output")
    assert not (tmp_path / "stale_output").exists()
    validate_cases({"cases": cases})
    with pytest.raises(ValueError):
        validate_cases({"cases": cases[:-1]})
    cases[0]["state"] = "legacy"
    with pytest.raises(ValueError):
        validate_cases({"cases": cases})
