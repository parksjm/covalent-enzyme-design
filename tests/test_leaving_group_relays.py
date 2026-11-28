"""Verify directed explicit-proton paths and independent topology mapping."""

from __future__ import annotations

from typing import Any

import numpy as np
import parmed
import pytest

from nyl12_covalent.leaving_group_relay_mapping import resolve_atom, resolve_roles
from nyl12_covalent.leaving_group_relays import (
    best_edge,
    catalytic_geometry,
    screen_paths,
)

CRITERIA = {"da_A": 3.5, "ha_A": 2.5, "angle_deg": 120}
CELL = np.eye(3) * 30


def role(label: str, index: int, hydrogen: int, water: bool = False) -> dict:
    """Return a small resolved synthetic chemical role."""
    return {
        "label": label,
        "index": index,
        "h_indices": [hydrogen],
        "water": water,
        "acceptor": water,
        "start": not water,
        "internal": water,
        "state": "H2O" if water else "OH",
        "conditional": not water,
    }


def test_two_water_route_and_periodic_translation() -> None:
    """A strict three-edge route survives wrapping; reverse H orientations fail."""
    xyz = np.array(
        [
            [0.0, 0, 0],
            [1.0, 0, 0],
            [3.0, 0, 0],
            [4.0, 0, 0],
            [6.0, 0, 0],
            [7.0, 0, 0],
            [9.0, 0, 0],
        ]
    )
    roles = [role("OG1", 0, 1), role("w1", 2, 3, True), role("w2", 4, 5, True)]
    rows = screen_paths(xyz, roles, 6, CELL, CRITERIA)
    assert len(rows) == 1 and rows[0]["strict"]
    assert rows[0]["water_count"] == 2
    assert rows[0]["path"] == ["OG1", "w1", "w2", "leaving_amide_N"]
    xyz[2:4] += CELL[0]
    assert screen_paths(xyz, roles, 6, CELL, CRITERIA)[0]["strict"]
    xyz[1] = [-1.0, 0, 0]
    failed = screen_paths(xyz, roles, 6, CELL, CRITERIA)
    assert not failed[0]["strict"]
    assert failed[0]["violation"] > 0


def test_no_hydrogen_donor_and_deprotonated_internal() -> None:
    """An anionic acceptor without an explicit proton cannot relay outward."""
    xyz = np.array([[0.0, 0, 0], [1.0, 0, 0], [3.0, 0, 0], [6.0, 0, 0]])
    asp = {**role("Asp", 2, 1, True), "h_indices": [], "water": False}
    rows = screen_paths(xyz, [role("OG1", 0, 1), asp], 3, CELL, CRITERIA)
    assert not rows[0]["strict"] and rows[0]["water_count"] == 0
    with pytest.raises(ValueError, match="no explicitly"):
        best_edge(xyz, asp, 3, CELL, CRITERIA)


def structure() -> Any:
    """Build one NH2 terminus and bonded water."""
    result = parmed.Structure()
    for name, z, residue, number in [
        ("N", 7, "THR", 1),
        ("H1", 1, "THR", 1),
        ("H2", 1, "THR", 1),
        ("O", 8, "WAT", 2),
        ("H1", 1, "WAT", 2),
        ("H2", 1, "WAT", 2),
    ]:
        result.add_atom(parmed.Atom(name=name, atomic_number=z), residue, number)
    for a, b in [(0, 1), (0, 2), (3, 4), (3, 5)]:
        result.bonds.append(parmed.Bond(result.atoms[a], result.atoms[b]))
    return result


def config() -> dict:
    """Return explicit mapping with neutral N-alpha state."""
    atom = {"residue_index1": 1, "resname": "THR", "atom": "N"}
    return {
        "roles": [
            {
                **atom,
                "label": "Nalpha",
                "hydrogens": ["H1", "H2"],
                "state": "NH2",
                "acceptor": True,
                "start": True,
                "internal": True,
                "conditional": True,
            }
        ],
        "geometry_atoms": {"nalpha": atom},
    }


def test_topology_mapping_and_charged_state_validation() -> None:
    """Explicit inventory and residue identity reject changed topology mistakes."""
    top, cfg = structure(), config()
    roles, atoms = resolve_roles(top, cfg)
    assert atoms["nalpha"] == 0 and len(roles) == 2
    assert roles[0]["h_indices"] == [1, 2]
    for update, message in [
        ({"resname": "LYS"}, "identity"),
        ({"chain": "Z"}, "Chain"),
        ({"atom": "NZ"}, "ambiguous"),
        ({"residue_index1": 0}, "outside"),
    ]:
        with pytest.raises(ValueError, match=message):
            resolve_atom(top, {**cfg["geometry_atoms"]["nalpha"], **update})
    cfg["roles"][0]["state"] = "NH3+"
    with pytest.raises(ValueError, match="NH3"):
        resolve_roles(top, cfg)
    top.add_atom(parmed.Atom(name="H3", atomic_number=1), "THR", 3)
    top.bonds.append(parmed.Bond(top.atoms[0], top.atoms[-1]))
    cfg["roles"][0].update(hydrogens=["H1", "H2", "H3"], acceptor=False)
    assert len(resolve_roles(top, cfg)[0][0]["h_indices"]) == 3
    cfg["roles"][0]["hydrogens"] = ["H1"]
    with pytest.raises(ValueError, match="inventory"):
        resolve_roles(top, cfg)


@pytest.mark.parametrize(
    "change,message",
    [
        ({"state": ""}, "unique"),
        ({"state": "COO-"}, "Deprotonated"),
        ({"state": "NH2", "hydrogens": []}, "inventory"),
    ],
)
def test_invalid_role_states(change: dict, message: str) -> None:
    """Invalid explicitly declared proton states fail rather than being repaired."""
    cfg = config()
    cfg["roles"][0].update(change)
    with pytest.raises(ValueError, match=message):
        resolve_roles(structure(), cfg)


def test_catalytic_geometry_is_fresh_and_periodic() -> None:
    """Reuse signed phi and catalytic distances after a common periodic unwrap."""
    xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -0.5, 0.0],
            [0.0, -1.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
        ]
    )
    mapped = dict(
        zip(
            [
                "carbon",
                "oxygen",
                "nleave",
                "adjacent_c",
                "nalpha",
                "og1",
                "oxyanion1",
                "oxyanion2",
            ],
            range(8),
            strict=True,
        )
    )
    before = catalytic_geometry(xyz, mapped, CELL)
    xyz[1] += CELL[0]
    after = catalytic_geometry(xyz, mapped, CELL)
    assert before == pytest.approx(after)
    assert before["d1"] == pytest.approx(np.sqrt(3))


def test_streaming_cli_and_refinement(tmp_path: Any, monkeypatch: Any) -> None:
    """CLI writes reusable route provenance and selects requested source frames."""
    import importlib.util
    import json
    from pathlib import Path
    import sys

    from scipy.io import netcdf_file

    script = (
        Path(__file__).resolve().parents[1] / "scripts/screen_leaving_group_relays.py"
    )
    spec = importlib.util.spec_from_file_location("leaving_cli", script)
    assert spec and spec.loader
    cli = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cli)
    top = structure()
    top.add_atom(parmed.Atom(name="N", atomic_number=7), "LIG", 3)
    topology = tmp_path / "top.prmtop"
    topology.write_text("Monkeypatched topology fixture")
    trajectory = tmp_path / "md.nc"
    xyz = np.array(
        [
            [0.0, 0, 0],
            [1.0, 0, 0],
            [-1.0, 0, 0],
            [3.0, 0, 0],
            [4.0, 0, 0],
            [3.0, 1, 0],
            [6.0, 0, 0],
        ]
    )
    with netcdf_file(str(trajectory), "w") as nc:
        for name, count in [("frame", 2), ("atom", 7), ("xyz", 3)]:
            nc.createDimension(name, count)
        nc.createVariable("coordinates", "f", ("frame", "atom", "xyz"))[:] = [xyz, xyz]
        nc.createVariable("cell_lengths", "f", ("frame", "xyz"))[:] = [[30] * 3] * 2
        nc.createVariable("cell_angles", "f", ("frame", "xyz"))[:] = [[90] * 3] * 2
        nc.createVariable("time", "f", ("frame",))[:] = [10, 20]
    case = {
        **config(),
        "name": "fixture",
        "topology": str(topology),
        "trajectories": [str(trajectory)],
    }
    case["geometry_atoms"]["nleave"] = {
        "residue_index1": 3,
        "resname": "LIG",
        "atom": "N",
    }
    monkeypatch.setattr(cli.parmed, "load_file", lambda _: top)
    monkeypatch.setattr(
        cli,
        "catalytic_geometry",
        lambda *_: {
            "d1": 3.0,
            "d2": 2.0,
            "d3": 2.0,
            "attack": 105.0,
            "projected_phi": -90.0,
        },
    )
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps({"stride": 5, "criteria": CRITERIA, "cases": [case]})
    )
    frames = tmp_path / "frames.json"
    frames.write_text(json.dumps([{"case": "fixture", "frame": 2}]))
    output = tmp_path / "screen"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(script),
            "--config",
            str(config_path),
            "--output",
            str(output),
            "--frames",
            str(frames),
        ],
    )
    cli.main()
    summary = json.loads((output / "summary.json").read_text())
    result = summary["cases"][0]
    assert result["counts"]["sampled_frames"] == 1
    assert result["water_count"] == 1
    assert len(result["role_mapping"]) == 1
    rows = [
        json.loads(line)
        for line in (output / "fixture_paths.jsonl").read_text().splitlines()
    ]
    assert rows[0]["frame"] == 2 and rows[0]["local_frame"] == 2
    assert rows[0]["projected_phi"] == -90
    with pytest.raises(ValueError, match="outside"):
        cli.scan_case(case, CRITERIA, tmp_path, 1, {3})
