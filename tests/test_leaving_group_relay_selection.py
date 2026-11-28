"""Selection gates, temporal diversity and source-coordinate preservation."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.io import netcdf_file

from nyl12_covalent import leaving_group_relay_export as export
from nyl12_covalent.leaving_group_relay_selection import (
    diverse,
    final_selection,
    read_candidates,
    refinement_frames,
    shortlist,
)


def row(frame: int = 1, case: str = "r1", **values: object) -> dict:
    result = dict(
        case=case,
        frame=frame,
        local_frame=frame,
        time_ps=frame * 10,
        start="OG",
        strict=True,
        d1=3.0,
        d2=3.0,
        d3=3.0,
        attack=105,
        projected_phi=-90,
        tier="A",
        native_pass=True,
        violation=0,
        water_count=1,
        internal_count=0,
        atom_indices1=[1, 4, 6],
        path=["Nalpha", "water:2", "leaving_amide_N"],
    )
    result.update(values)
    return result


def test_native_labels_cannot_change_relay_selection_order() -> None:
    """A diagnostic flag must not displace a chemically ranked candidate."""
    candidates = [row(1, native_pass=False), row(201, native_pass=True)]
    changed = [
        {**candidate, "native_pass": not candidate["native_pass"]}
        for candidate in candidates
    ]
    before = [candidate["frame"] for candidate in final_selection(candidates)]
    after = [candidate["frame"] for candidate in final_selection(changed)]
    assert before == after


def test_read_refinement_supersedes_and_preserves_phi(tmp_path: Path) -> None:
    directories = [tmp_path / "coarse", tmp_path / "fine"]
    for directory in directories:
        directory.mkdir()
    (directories[0] / "r_paths.jsonl").write_text(
        "\n".join(
            json.dumps(r)
            for r in [row(), row(2, strict=False), row(2, projected_phi=90)]
        )
        + "\n"
    )
    (directories[1] / "r_paths.jsonl").write_text(json.dumps(row(d1=8)) + "\n")
    role = dict(label="OG", state="OH", conditional=False, residue_index1=1, atom="OG1")
    config = {
        "cases": [
            {
                "name": "r1",
                "roles": [role],
                "geometry_atoms": {"og1": role, "nalpha": role},
            }
        ]
    }
    rows = read_candidates(directories, config)
    assert rows[0]["frame"] == 2
    assert rows[0]["phi_orientation"] == "positive"
    assert rows[1]["tier"] is None
    assert rows[1]["phi_orientation"] == "negative"
    with pytest.raises(ValueError, match="No .*screen records"):
        read_candidates([tmp_path / "absent"], config)


def test_diversity_and_bounded_selection() -> None:
    rows = [row(), row(2), row(101), row(102), row(201), row(301), row(case="r2")]
    selected = diverse(rows, 4)
    assert [(r["case"], r["frame"]) for r in selected] == [
        ("r1", 1),
        ("r2", 1),
        ("r1", 101),
        ("r1", 201),
    ]
    finals = final_selection(rows)
    assert [r["selection_role"] for r in finals] == [
        "primary",
        "primary",
        "reserve",
        "reserve",
    ]
    assert len(shortlist(rows, 3)) == 3
    assert final_selection([row(strict=False), row(tier=None)]) == []
    assert shortlist([row(strict=False)], 20)[0]["strict"] is False
    assert refinement_frames([row(), row(8)], {"r1": 10}) == [
        {"case": "r1", "frame": i} for i in range(1, 11)
    ]


def trajectory(path: Path) -> None:
    with netcdf_file(path, "w") as handle:
        for name, count in (("frame", 1), ("atom", 6), ("spatial", 3)):
            handle.createDimension(name, count)
        handle.createVariable("coordinates", "f", ("frame", "atom", "spatial"))[:] = (
            np.arange(18).reshape(1, 6, 3) / 7
        )
        handle.createVariable("time", "d", ("frame",))[:] = [10]
        for key, value in (
            ("cell_lengths", [30, 31, 32]),
            ("cell_angles", [90, 90, 90]),
        ):
            handle.createVariable(key, "d", ("frame", "spatial"))[:] = [value]


def topology() -> SimpleNamespace:
    residues = [
        SimpleNamespace(atoms=[SimpleNamespace(idx=i) for i in range(3)]),
        SimpleNamespace(atoms=[SimpleNamespace(idx=i) for i in range(3, 6)]),
    ]
    return SimpleNamespace(atoms=list(range(6)), residues=residues)


def test_export_exact_and_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, top = tmp_path / "source.nc", tmp_path / "top.prmtop"
    trajectory(path)
    top.write_text("immutable topology")
    candidate = row(trajectory=str(path))
    monkeypatch.setattr(export.parmed, "load_file", lambda _: topology())
    monkeypatch.setattr(export, "resolve_atom", lambda *_: 0)
    result = export.export_frame(
        candidate,
        {"topology": str(top), "geometry_atoms": {"nalpha": {}}},
        tmp_path,
        "primary",
    )
    assert result["coordinates_unmodified"]
    assert result["water_atoms_index1"] == {"2": [4, 5, 6]}
    assert result["geometry_atoms_index1"] == {"nalpha": 1}
    assert len(result["restart_sha256"]) == 64
    with pytest.raises(ValueError, match="time mismatch"):
        export.source_frame({**candidate, "time_ps": 500})
    monkeypatch.setattr(export.parmed, "load_file", lambda _: SimpleNamespace(atoms=[]))
    with pytest.raises(ValueError, match="atom count"):
        export.export_frame(candidate, {"topology": str(top)}, tmp_path, "bad")


def test_native_gate_failures_preserved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "source.nc"
    trajectory(path)
    config = {
        "native_thresholds": {
            "og1_shift_A": 0.75,
            "nalpha_shift_A": 1.25,
            "abs_psi_delta_deg": 35,
            "abs_chi1_delta_deg": 20,
        },
        "reference_cif": "ref.cif",
        "cases": [{"name": "r1", "topology": "top", "geometry_atoms": {"nalpha": {}}}],
    }
    monkeypatch.setattr(export, "_cif_atoms", lambda *_: {})
    monkeypatch.setattr(export.parmed, "load_file", lambda _: topology())
    monkeypatch.setattr(export, "resolve_atom", lambda *_: 0)
    metrics = dict(
        native_pass=False,
        native_og1_shift_A=1,
        native_nalpha_shift_A=2,
        psi_delta_deg=40,
        chi1_delta_deg=25,
    )
    monkeypatch.setattr(export, "native_geometry", lambda *_: metrics)
    checked = export.evaluate_native(
        [row(trajectory=str(path)), row(trajectory=str(path))], config
    )
    assert checked[0]["native_failed_gates"] == [
        "og1_shift",
        "nalpha_shift",
        "psi",
        "chi1",
    ]
    assert not final_selection(checked)[0]["native_pass"]
    config["native_thresholds"] = {}
    with pytest.raises(ValueError, match="thresholds differ"):
        export.evaluate_native([], config)


def test_context_ranking_and_activation_leg(monkeypatch: pytest.MonkeyPatch) -> None:
    from nyl12_covalent.leaving_group_relay_selection import annotate_context

    roles = [
        dict(label="OG", state="OH", conditional=False, residue_index1=1, atom="OG1"),
        dict(label="N", state="NH2", conditional=True, residue_index1=1, atom="N"),
        dict(
            label="Asp", state="COOH", conditional=False, residue_index1=2, atom="OD2"
        ),
        dict(label="Lys", state="NH3+", conditional=True, residue_index1=3, atom="NZ"),
    ]
    case = {
        "name": "r1",
        "topology": "top",
        "roles": roles,
        "geometry_atoms": {"og1": roles[0], "nalpha": roles[1]},
    }
    config = {
        "native_thresholds": {
            "og1_shift_A": 0.75,
            "nalpha_shift_A": 1.25,
            "abs_psi_delta_deg": 35,
            "abs_chi1_delta_deg": 20,
        },
        "reference_cif": "ref.cif",
        "cases": [case],
    }
    ranked = annotate_context(
        [row(start=label) for label in ("N", "OG", "Asp", "Lys")], config
    )
    assert [r["start"] for r in ranked] == ["OG", "Asp", "N", "Lys"]
    candidate = next(r for r in ranked if r["start"] == "N")
    assert candidate["nalpha_NH2_start"]
    xyz = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]])
    donor = SimpleNamespace(bond_partners=[SimpleNamespace(idx=1, atomic_number=1)])
    monkeypatch.setattr(
        export.parmed, "load_file", lambda _: SimpleNamespace(atoms=[donor])
    )
    monkeypatch.setattr(
        export,
        "resolve_atom",
        lambda _, descriptor: 0 if descriptor["atom"] == "OG1" else 2,
    )
    monkeypatch.setattr(
        export, "source_frame", lambda _: (xyz, np.array([30, 30, 30, 90, 90, 90]))
    )
    monkeypatch.setattr(export, "_cif_atoms", lambda *_: {})
    monkeypatch.setattr(
        export,
        "native_geometry",
        lambda *_: dict(
            native_pass=True,
            native_og1_shift_A=0,
            native_nalpha_shift_A=0,
            psi_delta_deg=0,
            chi1_delta_deg=0,
        ),
    )
    checked = export.evaluate_native([candidate], config)
    assert checked[0]["activation_status"] == "direct_activation_observed"
    xyz[2] = [5, 0, 0]
    assert (
        export.evaluate_native([candidate], config)[0]["activation_status"]
        == "activation_unresolved"
    )
    donor.bond_partners = []
    with pytest.raises(ValueError, match="explicit owned hydrogen"):
        export.evaluate_native([candidate], config)


def test_cli_diagnostic_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from scripts import select_leaving_group_relay_candidates as cli

    screen = tmp_path / "screen"
    screen.mkdir()
    role = dict(label="OG", state="OH", conditional=False, residue_index1=1, atom="OG1")
    config = {
        "cases": [
            {
                "name": "r1",
                "roles": [role],
                "geometry_atoms": {"og1": role, "nalpha": role},
            }
        ]
    }
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(config))
    (screen / "r_paths.jsonl").write_text(
        json.dumps(row(start="OG", strict=False)) + "\n"
    )
    (screen / "summary.json").write_text(
        json.dumps({"cases": [{"case": "r1", "total_frames": 100, "input_sha256": {}}]})
    )
    output = tmp_path / "selection"
    monkeypatch.setattr(cli, "evaluate_native", lambda rows, _: rows)
    monkeypatch.setattr(
        "sys.argv",
        [
            "select",
            "--config",
            str(cfg),
            "--screen-dir",
            str(screen),
            "--output",
            str(output),
        ],
    )
    cli.main()
    assert (
        json.loads((output / "selection.json").read_text())["status"]
        == "no_qualified_candidate"
    )
    assert len(json.loads((output / "refinement_frames.json").read_text())) == 5
    assert "diagnostic nearpasses only" in (output / "summary.md").read_text()


def test_route_context_precedes_frame_deduplication(tmp_path: Path) -> None:
    """A longer strict OG1 route survives a shorter strict conditional route."""
    og = dict(label="OG", state="OH", conditional=False, residue_index1=1, atom="OG1")
    nalpha = dict(label="N", state="NH2", conditional=True, residue_index1=1, atom="N")
    config = {
        "cases": [
            {
                "name": "r1",
                "roles": [og, nalpha],
                "geometry_atoms": {"og1": og, "nalpha": nalpha},
            }
        ]
    }
    shorter = row(start="N", water_count=0, atom_indices1=[1, 6])
    connected = row(start="OG", water_count=2, atom_indices1=[2, 3, 4, 6])
    source = tmp_path / "r_paths.jsonl"
    for order in ([shorter, connected], [connected, shorter]):
        source.write_text("".join(json.dumps(route) + "\n" for route in order))
        selected = read_candidates([tmp_path], config)
        assert len(selected) == 1
        assert selected[0]["start"] == "OG"
        assert selected[0]["network_context"] == "OG1_start_connected_path"
        assert not selected[0]["nalpha_NH2_start"]
    connected["strict"] = False
    source.write_text(
        "".join(json.dumps(route) + "\n" for route in [shorter, connected])
    )
    assert read_candidates([tmp_path], config)[0]["start"] == "N"
