"""Synthetic runtime and native artifact gates; no scientific engine execution."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import parmed
from parmed.amber import Rst7
import pytest

from nyl12_covalent import relay_case as case
from nyl12_covalent import relay_stage as stage
from nyl12_covalent.native_sander_cli import ANGSTROM_PER_BOHR, HARTREE_TO_KCAL_PER_MOL


@pytest.fixture
def explicit_case(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple:
    names = [
        "thr_n",
        "thr_h1",
        "thr_h2",
        "thr_cb",
        "thr_og",
        "proton",
        "carbon",
        "oxygen",
        "leaving_n",
        "leaving_h",
        "water_o",
        "water_h1",
        "water_h2",
        "ash_oxygen",
        "ash_proton",
        "asp_o1",
        "asp_o2",
        "phe_n",
        "asn_nd2",
    ]
    numbers = [7, 1, 1, 6, 8, 1, 6, 8, 7, 1, 8, 1, 1, 8, 1, 8, 8, 7, 7]
    numbers += [1] * (54 - numbers.count(1))
    numbers += [6] * (102 - len(numbers)) + [6, 6]
    topology = parmed.Structure()
    for index, number in enumerate(numbers):
        topology.add_atom(parmed.Atom(name=f"A{index}", atomic_number=number), "RES", 1)
    if (sum(numbers[:102]) + 1) % 2:
        topology.bonds.append(parmed.Bond(topology.atoms[0], topology.atoms[103]))
    restart = Rst7(natom=104)
    restart.coordinates = np.arange(312).reshape(104, 3) * 0.03
    restart.box = [40, 40, 40, 90, 90, 90]
    source = tmp_path / "source.ncrst"
    restart.write(str(source), netcdf=True)
    topology_path = tmp_path / "topology.prmtop"
    topology_path.write_text("synthetic topology handled by injected loader")
    calibration = tmp_path / "calibration.json"
    calibration.write_text("{}")
    amber = tmp_path / "amber"
    for name in ("bin/sander", "lib/libdlfind.so", "lib/libsander.so"):
        path = amber / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic binary identity")
    cfg = dict(
        topology=str(topology_path),
        source_restart=str(source),
        qm_atom_ids=list(range(1, 103)),
        active_atom_ids=list(range(1, 104)),
        atom_map={name: i + 1 for i, name in enumerate(names)},
        mdin="&cntrl maxcyc=1 /\n&qmmm qmmask='@1-102',qmcharge=-1,qm_theory='DFTB3',dftb_telec=0.0 /\n&dlfind active='@1-103',maxcycle=1 /\n",
        qm_charge=-1,
        calibration=str(calibration),
        amber_home=str(amber),
        max_atom_step_A=0.05,
        force_gate=0.12,
        memory=7,
    )
    path = tmp_path / "case.json"
    path.write_text(json.dumps(cfg))
    monkeypatch.setattr(case.parmed, "load_file", lambda _: topology)
    monkeypatch.setattr(case, "load_calibration", lambda _: (-2240.0, {}))
    return path, cfg, topology


def test_load_and_state(explicit_case: tuple) -> None:
    path, _, _ = explicit_case
    model = case.load_case(path)
    report = case.state_record(model, model.coordinates)
    assert len(model.qm) == 102 and len(report["proton_owners"]) == 54
    assert report["q"] == pytest.approx(report["CN_N"] - report["CN_OG1"])
    assert report["rCN"] > 0


@pytest.mark.parametrize(
    "change",
    [
        "selection",
        "mask",
        "onecycle",
        "count",
        "charge",
        "method",
        "map",
        "element",
        "cap",
        "memory",
        "parity",
    ],
)
def test_load_rejects(explicit_case: tuple, change: str) -> None:
    path, cfg, topology = explicit_case
    if change == "selection":
        cfg["qm_atom_ids"][0] = 0
    elif change == "mask":
        cfg["mdin"] = cfg["mdin"].replace("@1-102", "@2-102")
    elif change == "onecycle":
        cfg["mdin"] = cfg["mdin"].replace("maxcyc=1", "maxcyc=2")
    elif change == "count":
        cfg["qm_atom_ids"] = cfg["qm_atom_ids"][:-1]
    elif change == "charge":
        cfg["qm_charge"] = 0
        cfg["mdin"] = cfg["mdin"].replace("qmcharge=-1", "qmcharge=0")
    elif change == "method":
        cfg["mdin"] = cfg["mdin"].replace("DFTB3", "DFTB2")
    elif change == "map":
        cfg["atom_map"]["thr_n"] = 104
    elif change == "element":
        topology.atoms[0].atomic_number = 9  # preserve parity
    elif change == "cap":
        cfg["max_atom_step_A"] = 0
    elif change == "memory":
        cfg["memory"] = 0
    elif change == "parity":
        topology.atoms[100].atomic_number = 7
    path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError):
        case.load_case(path)


def test_capped_step_units() -> None:
    raw = np.array([3.0, 4.0, 0, 0, 0, 2])
    result = stage.capped_step(raw, 0.05).reshape(-1, 3)
    assert np.linalg.norm(result[0]) == pytest.approx(0.05)
    assert np.linalg.norm(result[1]) == pytest.approx(0.02)
    np.testing.assert_allclose(
        stage.capped_step(np.ones(3), 5), np.ones(3) * ANGSTROM_PER_BOHR
    )
    with pytest.raises(ValueError):
        stage.capped_step(np.array([np.nan, 0, 0]), 0.05)


@pytest.mark.parametrize("problem", [None, "scc", "energies", "geometry", "forces"])
def test_checked_native_outputs(tmp_path: Path, problem: str | None) -> None:
    class Native:
        counter = 1
        branch_dir = tmp_path
        evidence = [{"evaluation": 1}]

        def __call__(self, xyz: np.ndarray) -> tuple:
            return 1.0, np.zeros_like(xyz)

    xyz = np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    mdout = "Energy calculation finished, energy: 1.0\n"
    if problem == "scc":
        mdout += "SCC failed\n"
    if problem == "energies":
        mdout += mdout
    (tmp_path / "eval_00001.mdout").write_text(mdout)
    for suffix in ("path_force.xyz", "path_active.xyz"):
        points = xyz.copy()
        if problem == "geometry" and "active" in suffix:
            points[0, 0] += 0.01
        text = "2\nframe\n" + "".join(f"X {x} {y} {z}\n" for x, y, z in points)
        if problem == "forces" and "force" in suffix:
            text += text
        (tmp_path / f"eval_00001.{suffix}").write_text(text)
    physical = case.CheckedPhysical(Native(), np.ones(2, dtype=bool))
    if problem:
        with pytest.raises(ValueError):
            physical(xyz)
    else:
        energy, gradient = physical(xyz)
        assert energy == 1 and gradient.shape == xyz.shape
        assert (tmp_path / "native_evidence.jsonl").exists()


@pytest.mark.parametrize("mode", ["budget", "converged", "failure", "ascent"])
def test_synthetic_stage(
    explicit_case: tuple, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    path, cfg, _ = explicit_case
    calls = []

    class Native:
        def __init__(self, *args: object) -> None:
            self.counter = 0

        def __call__(self, xyz: np.ndarray) -> tuple:
            self.counter += 1
            if mode == "failure":
                raise ValueError("synthetic SCC failure")
            return 2.0, np.full_like(xyz, 0.01 if mode == "converged" else 1.0)

    class Optimizer:
        provenance = {"test": "synthetic native optimizer"}

        def __init__(self, *args: object, **kwargs: object) -> None:
            assert kwargs["memory"] == 7

        def __enter__(self) -> Optimizer:
            return self

        def __exit__(self, *args: object) -> None:
            pass

        def step(self, x: np.ndarray, g: np.ndarray) -> np.ndarray:
            calls.append((x.copy(), g.copy()))
            return g if mode == "ascent" else -g

        def restart(self) -> None:
            calls.append("reset")

    monkeypatch.setattr(stage, "NativeSanderEvaluator", Native)
    monkeypatch.setattr(stage, "CheckedPhysical", lambda native, _: native)
    monkeypatch.setattr(stage, "DLFindLBFGS", Optimizer)
    source = tmp_path / "lbfgs.F90"
    source.write_text("synthetic optimizer source")
    monkeypatch.setattr(stage, "LBF_SOURCE", source)
    output = tmp_path / "stage"
    if mode in ("failure", "ascent"):
        with pytest.raises(ValueError):
            stage.run_stage(path, Path(cfg["source_restart"]), output, [], max_steps=2)
        assert not json.loads((output / "process.json").read_text())["completed"]
        return
    result = stage.run_stage(path, Path(cfg["source_restart"]), output, [], max_steps=2)
    assert result["forceconv"] is (mode == "converged")
    assert Path(result["endpath"]).is_file()
    copied_names = [path.name for path in (output / "source_snapshot").iterdir()]
    for name in (
        "relay_stage.py",
        "relay_case.py",
        "relay_stage_units.py",
        "proton_coordination.py",
        "occupancy_bias.py",
        "prepare_dimer_continuation.py",
    ):
        assert sum(path.endswith("_" + name) for path in copied_names) == 1
    assert not any(
        path.endswith(("_sander", "_libsander.so", "_libdlfind.so"))
        for path in copied_names
    )
    manifest = json.loads((output / "source_manifest.json").read_text())
    assert len(manifest["binaries"]) == 3
    assert all(len(item["sha256"]) == 64 for item in manifest["binaries"])

    rows = [
        json.loads(line)
        for line in (output / "accepted_iterations.jsonl").read_text().splitlines()
    ]
    assert rows[0]["bias_energy_kcal_mol"] == 0
    assert rows[0]["physical_energy_kcal_mol"] == pytest.approx(
        2 * HARTREE_TO_KCAL_PER_MOL
    )
    if mode == "budget":
        assert any(isinstance(call, str) and call == "reset" for call in calls)
        np.testing.assert_allclose(
            calls[0][1], ANGSTROM_PER_BOHR / HARTREE_TO_KCAL_PER_MOL
        )
        final = Rst7.open(result["endpath"]).coordinates.reshape(-1, 3)
        original = Rst7.open(cfg["source_restart"]).coordinates.reshape(-1, 3)
        np.testing.assert_array_equal(final[103:], original[103:])
    assert (
        json.loads((output / "scientific_acceptance.json").read_text())[
            "chemical_acceptance"
        ]
        is None
    )
