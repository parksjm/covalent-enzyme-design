"""Direct controls never hide assisted hypotheses or become new relay progress."""

import json
from pathlib import Path

import pytest

from nyl12_covalent.leaving_group_relay_assisted import (
    CLASSES,
    assisted_class,
    audit_shortlist,
    read_assisted,
    select_assisted,
)


def configuration() -> dict:
    role = dict(label="OG", state="OH", conditional=False, residue_index1=1, atom="OG1")
    return {
        "cases": [
            {
                "name": "r1",
                "roles": [role],
                "geometry_atoms": {"og1": role, "nalpha": role},
            }
        ]
    }


def route(frame: int = 1, **changes: object) -> dict:
    value = dict(
        case="r1",
        frame=frame,
        local_frame=frame,
        time_ps=frame * 10,
        start="OG",
        strict=False,
        d1=3,
        d2=3,
        d3=3,
        attack=105,
        projected_phi=-90,
        violation=0.15,
        water_count=1,
        internal_count=0,
        atom_indices1=[1, 4, 6],
        native_pass=False,
    )
    value.update(changes)
    return value


def screen(directory: Path, rows: list[dict]) -> None:
    directory.mkdir()
    (directory / "r_paths.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in rows)
    )
    (directory / "summary.json").write_text(
        json.dumps(
            {
                "cases": [
                    {
                        "case": "r1",
                        "total_frames": 1000,
                        "input_sha256": {"trajectory": "digest"},
                    }
                ]
            }
        )
    )


def test_native_labels_cannot_change_assisted_selection(tmp_path: Path) -> None:
    """Native resemblance is independent of chemical qualification and rank."""
    directory = tmp_path / "source"
    screen(
        directory,
        [
            route(1, strict=True, native_pass=False),
            route(201, strict=True, native_pass=True),
        ],
    )
    rows, _, _ = read_assisted(
        [directory], configuration(), tmp_path / "retained.jsonl"
    )
    changed = [
        {**candidate, "native_pass": not candidate["native_pass"]} for candidate in rows
    ]
    before = [candidate["frame"] for candidate in select_assisted(rows)]
    after = [candidate["frame"] for candidate in select_assisted(changed)]
    assert before == after


def test_strict_direct_control_does_not_mask_assisted_nearpass(tmp_path: Path) -> None:
    directory = tmp_path / "coarse"
    screen(
        directory,
        [
            route(strict=True, water_count=0),
            route(),
            route(101, strict=True, d1=5),
            route(201, strict=True),
        ],
    )
    retained = tmp_path / "retained.jsonl"
    rows, provenance, totals = read_assisted([directory], configuration(), retained)
    assert [row["source_class"] for row in rows] == [CLASSES[0], CLASSES[1], CLASSES[2]]
    nearpass = next(row for row in rows if row["frame"] == 1)
    assert nearpass["water_count"] == 1 and not nearpass["strict"]
    assert nearpass["network_context"] == "OG1_start_connected_path"
    assert nearpass["phi_orientation"] == "negative"
    assert provenance["path_counts"]["direct_strict_tier_ab_paths"] == 1
    assert len(retained.read_text().splitlines()) == 4
    assert totals == {"r1": 1000}
    selected = select_assisted(rows)
    assert [row["selection_role"] for row in selected] == [
        "joint_qualified_assisted",
        "diagnostic_nearpass_not_qualified",
        "diagnostic_nearpass_not_qualified",
    ]


def test_refinement_supersedes_best_assisted_and_audit_is_bounded(
    tmp_path: Path,
) -> None:
    coarse, refined = tmp_path / "coarse", tmp_path / "refined"
    screen(coarse, [route()])
    screen(refined, [route(strict=True), route(strict=False, violation=0.05)])
    rows, _, _ = read_assisted(
        [coarse, refined], configuration(), tmp_path / "all.jsonl"
    )
    assert rows[0]["strict"]
    diagnostics = [route(i * 101, source_class=CLASSES[1]) for i in range(30)]
    diagnostics += [route(i * 101 + 1, source_class=CLASSES[2]) for i in range(30)]
    checked = audit_shortlist(diagnostics)
    assert len(checked) == 20
    assert {row["source_class"] for row in checked} == {CLASSES[1], CLASSES[2]}
    assert (
        assisted_class(route(tier=None, internal_count=1, water_count=0)) == CLASSES[3]
    )
    assert select_assisted([]) == []


def test_incomplete_screen_fails(tmp_path: Path) -> None:
    source = tmp_path / "screen"
    source.mkdir()
    (source / "r_paths.jsonl").write_text(json.dumps(route()) + "\n")
    with pytest.raises(ValueError, match="Incomplete screen"):
        read_assisted([source], configuration(), tmp_path / "retained.jsonl")
    with pytest.raises(ValueError, match="No complete"):
        read_assisted([tmp_path / "absent"], configuration(), tmp_path / "empty.jsonl")


def test_cli_labels_diagnostics_and_writes_refinement_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from scripts import review_assisted_leaving_group_relays as cli

    source, config, output = (
        tmp_path / "coarse",
        tmp_path / "config.json",
        tmp_path / "review",
    )
    screen(source, [route(strict=True, water_count=0), route()])
    config.write_text(json.dumps(configuration()))
    monkeypatch.setattr(
        cli,
        "evaluate_native",
        lambda rows, _: [
            {**row, "native_failed_gates": ["chi1"], "native_pass": False}
            for row in rows
        ],
    )
    monkeypatch.setattr(cli, "export_frame", lambda row, *_: row)
    monkeypatch.setattr(
        "sys.argv",
        [
            "review",
            "--config",
            str(config),
            "--screen-dir",
            str(source),
            "--output",
            str(output),
        ],
    )
    cli.main()
    report = json.loads((output / "selection.json").read_text())
    assert report["status"] == "no_joint_qualified_assisted_candidate"
    assert report["exports"][0]["selection_role"] == "diagnostic_nearpass_not_qualified"
    assert report["exports"][0]["native_failed_gates"] == ["chi1"]
    assert not report["direct_controls_are_new_relays"]
    assert len(json.loads((output / "refinement_frames.json").read_text())) == 5
    assert json.loads((output / "refinement_cases.json").read_text())["stride"] == 1
