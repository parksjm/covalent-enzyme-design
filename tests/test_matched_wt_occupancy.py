"""Synthetic scientific and streaming checks for matched WT occupancy."""

from collections import Counter
import copy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.io import netcdf_file

from nyl12_covalent import matched_wt_occupancy as stats
from nyl12_covalent import matched_wt_occupancy_io as io
from nyl12_covalent.leaving_group_relays import catalytic_geometry


def config() -> dict:
    """Return a complete synthetic matched-cohort configuration."""
    keys = (
        "og1",
        "nalpha",
        "carbon",
        "oxygen",
        "nleave",
        "adjacent_c",
        "oxyanion1",
        "oxyanion2",
    )
    return {
        "dt_ps": 10,
        "time_offset_ps": 500,
        "native_gates": stats.NATIVE_GATES,
        "cases": [
            {
                "name": f"{state}_{i}",
                "state": state,
                "expected_frames": 10000,
                "trajectories": [f"{state}_{i}.nc"],
                "geometry_atoms": dict.fromkeys(keys, {}),
            }
            for state in sorted(stats.STATES)
            for i in range(4)
        ],
    }


def test_config_requires_matched_independent_complete_cohort() -> None:
    stats.validate_config(config())
    for mutate in (
        lambda c: c["cases"].pop(),
        lambda c: c.update(dt_ps=20),
        lambda c: c.update(time_offset_ps=0),
        lambda c: c["cases"][0].update(name="../escape"),
        lambda c: c["cases"][0].update(expected_frames=9999),
        lambda c: c["cases"][0].update(geometry_atoms={}),
        lambda c: c["cases"][1].update(trajectories=c["cases"][0]["trajectories"]),
    ):
        candidate = copy.deepcopy(config())
        mutate(candidate)
        with pytest.raises(ValueError):
            stats.validate_config(candidate)


def test_thresholds_exclusive_tiers_and_independent_phi() -> None:
    metrics = dict(d1=3.2, d2=3.8, d3=3.8, attack=90.0, projected_phi=70.0)
    positive = stats.frame_flags(metrics, True)
    assert positive["tier_a"] and not positive["tier_b_exclusive"]
    assert positive["tier_a_or_b_and_native"]
    negative = stats.frame_flags({**metrics, "projected_phi": -70.0}, False)
    assert negative["tier_a"] and negative["phi_negative"]
    assert not negative["tier_a_phi_positive"]
    relaxed = stats.frame_flags({**metrics, "d1": 3.5}, True)
    assert relaxed["tier_b_exclusive"] and relaxed["tier_a_or_b"]
    assert not relaxed["tier_a"]
    outside = stats.frame_flags({**metrics, "d1": 3.50001, "projected_phi": 0}, True)
    assert outside["phi_zero"] and not outside["tier_a_or_b"]
    with pytest.raises(ValueError, match="Nonfinite"):
        stats.frame_flags({**metrics, "d1": float("nan")}, False)


def test_production_block_endpoints_and_replica_not_frame_statistics() -> None:
    assert [stats.block_index(f, 10) for f in (1, 1000, 1001, 10000)] == [1, 1, 2, 10]
    with pytest.raises(ValueError):
        stats.block_index(0, 10)
    with pytest.raises(ValueError):
        stats.fractions(Counter(), 0)
    replicas = [
        {"state": state, "fractions": dict.fromkeys(stats.OBSERVABLES, value)}
        for state in stats.STATES
        for value in (0.0, 0.0, 1.0, 1.0)
    ]
    result = stats.summarize_states(replicas)
    row = result[next(iter(stats.STATES))]["observables"]["tier_a"]
    assert row == {
        "mean": 0.5,
        "sample_sd": pytest.approx(np.sqrt(1 / 3)),
        "min": 0.0,
        "max": 1.0,
    }
    with pytest.raises(ValueError):
        stats.summarize_states(replicas[:-1])
    assert stats.fractions(Counter(tier_a=3), 4)["tier_a"] == 0.75


def topology() -> tuple[SimpleNamespace, dict]:
    """Return an explicit independently indexed bonded miniature topology."""
    keys = (
        "og1",
        "nalpha",
        "carbon",
        "oxygen",
        "nleave",
        "adjacent_c",
        "oxyanion1",
        "oxyanion2",
    )
    atoms = [
        SimpleNamespace(idx=i, name=k, bond_partners=[]) for i, k in enumerate(keys)
    ]
    atoms[2].bond_partners = [atoms[i] for i in (3, 4, 5)]
    top = SimpleNamespace(
        atoms=atoms, residues=[SimpleNamespace(name="SYN", atoms=atoms)]
    )
    descriptors = {
        key: {"residue_index1": 1, "resname": "SYN", "atom": key} for key in keys
    }
    return top, descriptors


def test_mapping_connectivity_identity_and_pbc() -> None:
    top, descriptors = topology()
    mapped = io.resolve_geometry(top, descriptors)
    xyz = np.array(
        [
            [0, 2, 2],
            [0, 1, 0],
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [-1, -1, 0],
            [1, 3, 0],
            [1, 0, 3],
        ],
        dtype=float,
    )
    cell = np.eye(3) * 20
    first = catalytic_geometry(xyz, mapped, cell)
    xyz[[0, 3, 4]] += np.array([20, -20, 20])
    assert catalytic_geometry(xyz, mapped, cell) == pytest.approx(first)
    top.atoms[2].bond_partners.pop()
    with pytest.raises(ValueError, match="bond"):
        io.resolve_geometry(top, descriptors)
    descriptors["og1"]["resname"] = "WRONG"
    with pytest.raises(ValueError, match="identity"):
        io.resolve_geometry(top, descriptors)


def trajectory(path: Path, times: list[float]) -> None:
    """Write a tiny real Amber-like NetCDF stream for read-only reader checks."""
    with netcdf_file(str(path), "w") as source:
        source.createDimension("frame", len(times))
        source.createDimension("atom", 8)
        source.createDimension("spatial", 3)
        source.createVariable("coordinates", "f", ("frame", "atom", "spatial"))[:] = (
            np.zeros((len(times), 8, 3))
        )
        source.createVariable("cell_lengths", "f", ("frame", "spatial"))[:] = 20
        source.createVariable("cell_angles", "f", ("frame", "spatial"))[:] = 90
        source.createVariable("time", "d", ("frame",))[:] = times


def test_stream_all_frames_native_csv_hashes_and_failure_guards(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    top, descriptors = topology()
    top_path, ref = tmp_path / "top", tmp_path / "ref"
    top_path.write_text("topology")
    ref.write_text("reference")
    nc1, nc2 = tmp_path / "first.nc", tmp_path / "second.nc"
    trajectory(nc1, [510, 520])
    trajectory(nc2, [530, 540])
    monkeypatch.setattr(io.parmed, "load_file", lambda _: top)
    monkeypatch.setattr(io, "verify_microstate", lambda *_: None)
    monkeypatch.setattr(io, "_cif_atoms", lambda *args: {})
    metrics = dict(d1=3.0, d2=3.0, d3=3.0, attack=100.0, projected_phi=-90.0)
    monkeypatch.setattr(io, "catalytic_geometry", lambda *args: metrics)
    monkeypatch.setattr(
        io,
        "native_geometry",
        lambda *args: {"native_pass": True, "native_og1_shift_A": 0.1},
    )
    case = {
        "name": "replica",
        "state": "synthetic",
        "topology": str(top_path),
        "trajectories": [str(nc1), str(nc2)],
        "expected_frames": 4,
        "geometry_atoms": descriptors,
    }
    cfg = {"reference_cif": str(ref), "dt_ps": 10, "time_offset_ps": 500}
    out = tmp_path / "out"
    out.mkdir()
    result = io.analyze_case(case, cfg, out)
    assert result["frames"] == 4 and result["fractions"]["tier_a_and_native"] == 1
    assert result["counts"]["phi_negative"] == 4
    assert result["blocks_10ns"]["1"]["frames"] == 4
    assert result["input_sha256"][str(top_path)] == io.digest(top_path)
    assert len((out / "replica.csv").read_text().splitlines()) == 5
    for frames, offset, pattern in (
        (3, 500, "Too many"),
        (5, 500, "Missing"),
        (4, 0, "cadence"),
    ):
        attempt = tmp_path / f"bad{frames}_{offset}"
        attempt.mkdir()
        with pytest.raises(ValueError, match=pattern):
            io.analyze_case(
                {**case, "expected_frames": frames},
                {**cfg, "time_offset_ps": offset},
                attempt,
            )

    def change_source(*args: object) -> dict:
        ref.write_text(ref.read_text() + "changed")
        return {"native_pass": True}

    monkeypatch.setattr(io, "native_geometry", change_source)
    changed = tmp_path / "changed"
    changed.mkdir()
    with pytest.raises(ValueError, match="Source changed"):
        io.analyze_case(case, cfg, changed)


def test_cli_writes_new_summary_and_plot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import json
    import sys
    from scripts import analyze_matched_wt_occupancy as cli

    cfg = config()
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(json.dumps(cfg))
    output = tmp_path / "analysis"
    fractions = stats.fractions(Counter(tier_a=1, tier_a_or_b=1), 2)
    monkeypatch.setattr(
        cli,
        "analyze_case",
        lambda case, *_: {
            "name": case["name"],
            "state": case["state"],
            "fractions": fractions,
        },
    )
    monkeypatch.setattr(
        sys, "argv", ["analyze", "--config", str(cfg_path), "--output", str(output)]
    )
    cli.main()
    summary = json.loads((output / "summary.json").read_text())
    assert len(summary["replicas"]) == 12
    assert len(summary["source_sha256"]) == 10
    assert (output / "occupancy.png").stat().st_size > 1000
    with pytest.raises(FileExistsError):
        cli.main()


def test_declared_proton_inventory_all_four_sites() -> None:
    residues = [SimpleNamespace(name="unused", atoms=[]) for _ in range(1258)]
    atoms = []
    for thr, a, b in (
        (929, 968, 970),
        (1025, 1064, 1066),
        (1121, 1160, 1162),
        (1217, 1256, 1258),
    ):
        n = SimpleNamespace(
            idx=len(atoms),
            name="N",
            bond_partners=[SimpleNamespace(atomic_number=1)] * 2,
        )
        atoms.append(n)
        residues[thr - 1] = SimpleNamespace(name="THR", atoms=[n])
        for residue in (a, b):
            residues[residue - 1] = SimpleNamespace(name="ASP", atoms=[])
    top = SimpleNamespace(atoms=atoms, residues=residues)
    io.verify_microstate(top, "nh2_asp272_asp274")
    with pytest.raises(ValueError, match="Thr terminal"):
        io.verify_microstate(top, "nh3_asp272_asp274")
    with pytest.raises(ValueError, match="Asp residue"):
        io.verify_microstate(top, "nh2_asp272_ash274")
    for b in (970, 1066, 1162, 1258):
        residues[b - 1] = SimpleNamespace(
            name="ASH", atoms=[SimpleNamespace(name="HD2", atomic_number=1)]
        )
    io.verify_microstate(top, "nh2_asp272_ash274")
    residues[1258 - 1].atoms = []
    with pytest.raises(ValueError, match="Asp carboxyl"):
        io.verify_microstate(top, "nh2_asp272_ash274")
