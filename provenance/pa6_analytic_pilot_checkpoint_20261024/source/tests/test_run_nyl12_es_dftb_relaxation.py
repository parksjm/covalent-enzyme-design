"""Serialization, force-proof gating and mocked process tests; no engine runs."""

from __future__ import annotations

import copy
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_array_equal
import parmed
import pytest
from parmed.amber.mask import AmberMask
from parmed.topologyobjects import Atom, AtomType

from nyl12_covalent import es_dftb_relaxation as module


@pytest.mark.parametrize(
    ("indices", "expected"),
    [([1], "@1"), ([8, 1, 2, 7, 4], "@1-2,4,7-8")],
)
def test_atom_mask_ranges_preserve_gaps(indices: list[int], expected: str) -> None:
    """Range compression must not select atoms in an omitted gap."""
    assert module._compact_mask(indices) == expected


@pytest.mark.parametrize("indices", [[], [0], [True], [1, 1], list(range(1, 500, 2))])
def test_invalid_or_overlong_atom_mask_rejected(indices: list[int]) -> None:
    """Reject ambiguous atoms and masks that Amber would truncate."""
    with pytest.raises(ValueError):
        module._compact_mask(indices)


def _write(path: Path, value: dict | list) -> None:
    path.write_text(json.dumps(value))


@pytest.fixture
def prepared_inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict:
    """Test exact restart/restraint handling with an independently mocked core audit.

    Chemical validation is covered by test_es_qm_region; this fixture exercises
    the driver interface using98selectedatoms plus8unselectedhydrogens.
    """
    structure = parmed.Structure()
    for index in range(106):
        z = 6 if index < 98 else 1
        mass = 12.011 if z == 6 else 1.008
        typ = AtomType(f"T{z}", z, mass, z)
        typ.set_lj_params(0.1, 1.5)
        atom = Atom(name=f"A{index}", type=f"T{z}", atomic_number=z, mass=mass)
        atom.atom_type = typ
        structure.add_atom(atom, "TST", 1)
    structure.coordinates = np.random.default_rng(8).normal(size=(106, 3))
    structure.box = [81.123456, 89.234567, 93.345678, 90, 90, 90]
    parm = parmed.amber.AmberParm.from_structure(structure)
    topology, coordinates = tmp_path / "system.prmtop", tmp_path / "frame.ncrst"
    parm.save(str(topology))
    parm.save(str(coordinates), format="ncrst")
    amberhome = tmp_path / "amber"
    (amberhome / "bin").mkdir(parents=True)
    sander = amberhome / "bin/sander"
    sander.write_text("not executed")
    sander.chmod(0o755)
    python = tmp_path / "python"
    python.write_text("pinned Python executable")
    parameter = tmp_path / "C-H.skf"
    parameter.write_text("pinned parameter")
    api = tmp_path / "api.so"
    api.write_text("pinned compiled API identity")
    audit = {
        "qmmask": "@" + ",".join(map(str, range(1, 99))),
        "real_atom_count": 98,
        "link_H_count": 8,
        "electron_count_with_links": 368,
        "qm_atoms": [{"index": i} for i in range(1, 99)],
        "boundaries": [{"qm_index": i, "mm_index": 98 + i} for i in range(1, 9)],
    }
    monkeypatch.setattr(
        module, "audit_es_qm_region", lambda *args: copy.deepcopy(audit)
    )
    mm = {
        "ntb": 1,
        "igb": 0,
        "ntc": 1,
        "ntf": 1,
        "ntr": 0,
        "ifqnt": 1,
        "cut": 10.0,
        "jfastw": 4,
    }
    qm = {
        "qm_theory": "DFTB3",
        "qmcharge": -2,
        "qmmm_int": 1,
        "adjust_q": 2,
        "qm_ewald": 1,
        "qm_pme": 1,
        "qmcut": 10.0,
        "scfconv": 1e-10,
        "lnk_dis": 1.09,
        "lnk_atomic_no": 1,
        "lnk_method": 1,
        "qmmask": audit["qmmask"],
    }
    inventory = {
        "real_atoms": 98,
        "link_atoms": 8,
        "electron_count": 368,
        "boundary_pairs": [[i, 98 + i] for i in range(1, 9)],
    }
    engine = {
        "topology": str(topology),
        "coordinates": str(coordinates),
        "python_executable": str(python),
        "parameter_files": [str(parameter)],
        "environment": {
            "AMBERHOME": str(amberhome),
            "LD_LIBRARY_PATH": str(amberhome / "lib"),
            "OMP_NUM_THREADS": "1",
        },
        "mm_options": mm,
        "qm_options": qm,
        "expected_inventory": inventory,
        "components": [
            {"atom": i + 1, "axis": 0, "label": label}
            for i, label in enumerate(
                [
                    "Thr_OG1",
                    "PA6_C54",
                    "PA6_O49",
                    "QM_link_parent",
                    "MM_link_partner",
                    "external_MM_water_O",
                ]
            )
        ],
        "force_absolute_tolerance": 0.001,
        "force_relative_tolerance": 0.01,
        "step_sizes_angstrom": [0.0001, 0.0002],
        "energy_reproducibility_tolerance": 1e-6,
    }
    engine_path, region_path = tmp_path / "engine.json", tmp_path / "region.json"
    _write(engine_path, engine)
    _write(
        region_path,
        {
            "topology": str(topology),
            "coordinates": str(coordinates),
            "variants": [{"name": "primary_106_nuclei"}],
        },
    )
    engine["approved_inventory_source"] = {
        "path": str(region_path),
        "sha256": module._sha(region_path),
    }
    _write(engine_path, engine)
    force_dir = tmp_path / "force"
    (force_dir / "source").mkdir(parents=True)
    shutil.copy2(engine_path, force_dir / "source/declared_config.json")
    records = [
        {"path": str(path), "sha256": module._sha(path), "bytes": path.stat().st_size}
        for path in [engine_path, topology, coordinates, python, parameter]
    ]
    _write(force_dir / "source_manifest.json", records)
    full_forces = np.zeros((106, 3))
    full_forces[:6, 0] = 1.0
    baseline = {
        "inventory": inventory,
        "realized_options": {"qm": qm, "mm": mm},
        "energy_kcal_per_mol": 100.0,
        "box": parm.box.tolist(),
        "perturbation": None,
        "forces": [1.0] * 6,
        "api_sources": [
            {"path": str(api), "sha256": module._sha(api), "bytes": api.stat().st_size}
        ],
    }

    def worker(name: str, perturbation: dict | None, energy: float) -> None:
        directory = force_dir / name
        directory.mkdir()
        _write(
            directory / "request.json", {"config": engine, "perturbation": perturbation}
        )
        _write(
            directory / "result.json",
            {**baseline, "perturbation": perturbation, "energy_kcal_per_mol": energy},
        )
        _write(directory / "process_status.json", {"exit_status": 0})

    worker("baseline", None, 100.0)
    worker("baseline_repeat", None, 100.0 + 1e-8)
    np.save(force_dir / "baseline/analytic_forces.npy", full_forces)
    comparisons = []
    for i, probe in enumerate(engine["components"]):
        for j, step in enumerate(engine["step_sizes_angstrom"]):
            energies = [100.0 - sign * step * 1.001 for sign in [1, -1]]
            for name, sign, energy in zip(
                ["plus", "minus"], [1, -1], energies, strict=True
            ):
                worker(
                    f"component_{i}_step_{j}_{name}",
                    {
                        "atom": probe["atom"],
                        "axis": probe["axis"],
                        "delta": sign * step,
                    },
                    energy,
                )
            fd = -(energies[0] - energies[1]) / (2 * step)
            comparisons.append(
                {
                    **probe,
                    "step_angstrom": step,
                    "passed": True,
                    "analytic_force": 1.0,
                    "finite_difference_force": fd,
                    "absolute_error": abs(1.0 - fd),
                    "tolerance": 0.001 + 0.01 * max(1.0, abs(fd)),
                }
            )
    check = {
        "accepted": True,
        "sources_unchanged": True,
        "energy_repeat_difference": abs(100.0 - (100.0 + 1e-8)),
        "baseline": baseline,
        "components": comparisons,
    }
    force_path = force_dir / "force_check.json"
    _write(force_path, check)
    analytic = tmp_path / "analytic"
    analytic.mkdir()
    for name in ["result.json", "request.json", "analytic_forces.npy"]:
        shutil.copy2(force_dir / "baseline" / name, analytic / name)
    decision_path = tmp_path / "decision.json"
    decision = {
        "mode": "bounded_analytic_gradient_pilot",
        "numerical_force_check_required": False,
        "rationale": "User steered a bounded analytic-gradient pilot; numerical check remains failed.",
        "user_steering_reference": "Explicit test authorization",
        "baseline_directory": str(analytic),
        "baseline_result_sha256": module._sha(analytic / "result.json"),
        "baseline_forces_sha256": module._sha(analytic / "analytic_forces.npy"),
        "retained_failed_force_check": {
            "path": str(force_path),
            "sha256": module._sha(force_path),
        },
        "approved_inventory_source": engine["approved_inventory_source"],
        "max_cycles": 25,
        "mm_restraint_weight": 10.0,
        "reaction_coordinate_restraints": False,
        "scope": "Exploratory ES relaxation only; no chemical-model validation",
    }
    _write(decision_path, decision)
    return {
        "decision": decision_path,
        "engine": engine_path,
        "force": force_path,
        "region": region_path,
        "sander": sander,
        "output": tmp_path / "output",
        "parameter": parameter,
    }


def _prepare(inputs: dict, pilot: bool = False) -> dict:
    return module.prepare_es_relaxation(
        inputs["engine"],
        inputs["force"],
        inputs["region"],
        inputs["sander"],
        inputs["output"],
        inputs["decision"] if pilot else None,
    )


def test_exact_restart_all_mm_restraints_no_engine(
    prepared_inputs: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        module.subprocess,
        "run",
        lambda *args, **kwargs: pytest.fail("Preparation must not run engine"),
    )
    data = _prepare(prepared_inputs)
    output = prepared_inputs["output"]
    original = json.loads(prepared_inputs["engine"].read_text())
    before = parmed.amber.AmberParm(original["topology"], xyz=original["coordinates"])
    after = parmed.amber.AmberParm(
        original["topology"], xyz=str(output / "starting.ncrst")
    )
    assert_array_equal(before.coordinates, after.coordinates)
    assert_array_equal(before.box, after.box)
    mask = "!(@1-98)"
    selection = AmberMask(after, mask).Selection()
    assert selection[:98] == [0] * 98 and selection[98:] == [1] * 8
    mdin = (output / "relaxation.mdin").read_text()
    assert "maxcyc=25" in mdin and "restraint_wt=10.0" in mdin and "ntwx=1" in mdin
    assert "-ref" in data["command"] and "optimization.nc" in " ".join(data["command"])
    assert data["restraint"]["atom_count"] == 8
    assert not (output / "process.json").exists()
    with pytest.raises(FileExistsError):
        _prepare(prepared_inputs)


@pytest.mark.parametrize(
    "mutation",
    [
        "accepted",
        "sources",
        "component",
        "numeric",
        "inventory",
        "realized",
        "parameter",
        "snapshot",
        "physics",
    ],
)
def test_force_gate_rejects_invalid_evidence(
    prepared_inputs: dict, mutation: str
) -> None:
    path = prepared_inputs["force"]
    check = json.loads(path.read_text())
    if mutation == "accepted":
        check["accepted"] = False
    elif mutation == "sources":
        check["sources_unchanged"] = False
    elif mutation == "component":
        check["components"][0]["passed"] = False
    elif mutation == "numeric":
        check["components"][0]["analytic_force"] = 100
    elif mutation == "inventory":
        check["baseline"]["inventory"]["real_atoms"] = 97
    elif mutation == "realized":
        check["baseline"]["realized_options"]["qm"]["qm_pme"] = 0
    elif mutation == "parameter":
        prepared_inputs["parameter"].write_text("changed")
    elif mutation == "snapshot":
        (path.parent / "source/declared_config.json").write_text("{}")
    else:
        engine = json.loads(prepared_inputs["engine"].read_text())
        engine["qm_options"]["adjust_q"] = 1
        _write(prepared_inputs["engine"], engine)
        shutil.copy2(
            prepared_inputs["engine"], path.parent / "source/declared_config.json"
        )
        records = json.loads((path.parent / "source_manifest.json").read_text())
        for row in records:
            if row["path"] == str(prepared_inputs["engine"]):
                row.update(
                    sha256=module._sha(prepared_inputs["engine"]),
                    bytes=prepared_inputs["engine"].stat().st_size,
                )
        _write(path.parent / "source_manifest.json", records)
    _write(path, check)
    with pytest.raises(ValueError):
        _prepare(prepared_inputs)
    assert not prepared_inputs["output"].exists()


@pytest.mark.parametrize("exit_status", [0, 2])
def test_mocked_execution_records_separate_acceptance(
    prepared_inputs: dict, monkeypatch: pytest.MonkeyPatch, exit_status: int
) -> None:
    data = _prepare(prepared_inputs)
    output = prepared_inputs["output"]

    def fake(command: list, **kwargs: object) -> SimpleNamespace:
        shutil.copy2(output / "starting.ncrst", output / "optimized.ncrst")
        return SimpleNamespace(returncode=exit_status)

    monkeypatch.setattr(module.subprocess, "run", fake)
    if exit_status:
        with pytest.raises(RuntimeError):
            module.run_prepared_relaxation(data)
    else:
        assert module.run_prepared_relaxation(data)["exit_status"] == 0
    assert (
        json.loads((output / "process.json").read_text())["exit_status"] == exit_status
    )
    assert (
        json.loads((output / "geometry_acceptance.json").read_text())["accepted"]
        is False
    )
    with pytest.raises(ValueError):
        module.run_prepared_relaxation(data)


def test_post_prepare_parameter_mutation_blocks_execution(
    prepared_inputs: dict,
) -> None:
    data = _prepare(prepared_inputs)
    prepared_inputs["parameter"].write_text("changed after preparation")
    with pytest.raises(ValueError, match="source changed"):
        module.run_prepared_relaxation(data)


def _failed_decision(inputs: dict) -> dict:
    check = json.loads(inputs["force"].read_text())
    check["accepted"] = False
    _write(inputs["force"], check)
    decision = json.loads(inputs["decision"].read_text())
    decision["retained_failed_force_check"]["sha256"] = module._sha(inputs["force"])
    _write(inputs["decision"], decision)
    return decision


def test_explicit_analytic_pilot_preserves_failed_fd(prepared_inputs: dict) -> None:
    _failed_decision(prepared_inputs)
    data = _prepare(prepared_inputs, pilot=True)
    assert data["pilot_decision"]["mode"] == "bounded_analytic_gradient_pilot"
    assert data["numerical_force_validation_passed"] is False
    assert json.loads(prepared_inputs["force"].read_text())["accepted"] is False
    assert not (prepared_inputs["output"] / "process.json").exists()
    assert data["analytic_baseline_energy_kcal_per_mol"] == 100.0
    for filename in [
        "es_qm_region.py",
        "es_dftb_relaxation_gate.py",
        "es_dftb_analytic_pilot.py",
        "es_dftb_force_evidence.py",
    ]:
        assert (prepared_inputs["output"] / "source" / filename).is_file()


@pytest.mark.parametrize(
    "mutation",
    [
        "mode",
        "scope",
        "cycles",
        "restraint",
        "rationale",
        "fd_hash",
        "inventory",
        "baseline_hash",
        "request",
        "nonfinite",
        "realized_mm",
        "realized_qm",
    ],
)
def test_conditional_pilot_rejects_invalid_decisions(
    prepared_inputs: dict, mutation: str
) -> None:
    decision = _failed_decision(prepared_inputs)
    baseline = Path(decision["baseline_directory"])
    if mutation in {"mode", "scope", "rationale"}:
        decision[
            {"mode": "mode", "scope": "scope", "rationale": "rationale"}[mutation]
        ] = ""
    elif mutation == "cycles":
        decision["max_cycles"] = 100
    elif mutation == "restraint":
        decision["mm_restraint_weight"] = 0
    elif mutation == "fd_hash":
        decision["retained_failed_force_check"]["sha256"] = "bad"
    elif mutation == "inventory":
        decision["approved_inventory_source"]["sha256"] = "bad"
    elif mutation == "baseline_hash":
        decision["baseline_result_sha256"] = "bad"
    elif mutation == "request":
        request = json.loads((baseline / "request.json").read_text())
        request["config"]["qm_options"]["qmcharge"] = 0
        _write(baseline / "request.json", request)
    elif mutation == "nonfinite":
        forces = np.load(baseline / "analytic_forces.npy")
        forces[0, 0] = np.nan
        np.save(baseline / "analytic_forces.npy", forces)
        decision["baseline_forces_sha256"] = module._sha(
            baseline / "analytic_forces.npy"
        )
    else:
        result = json.loads((baseline / "result.json").read_text())
        category, option = (
            ("mm", "ntc") if mutation == "realized_mm" else ("qm", "qmcharge")
        )
        result["realized_options"][category][option] = 2 if category == "mm" else 0
        _write(baseline / "result.json", result)
        decision["baseline_result_sha256"] = module._sha(baseline / "result.json")
    _write(prepared_inputs["decision"], decision)
    with pytest.raises(ValueError):
        _prepare(prepared_inputs, pilot=True)
    assert not prepared_inputs["output"].exists()


@pytest.mark.parametrize(
    "mutation",
    ["omitted", "duplicate", "fabricated_tolerance", "worker_energy", "worker_status"],
)
def test_strict_gate_recomputes_preserved_proof(
    prepared_inputs: dict, mutation: str
) -> None:
    path = prepared_inputs["force"]
    check = json.loads(path.read_text())
    if mutation == "omitted":
        check["components"].pop()
    elif mutation == "duplicate":
        check["components"][-1] = check["components"][0]
    elif mutation == "fabricated_tolerance":
        check["components"][0]["tolerance"] = 100
    elif mutation == "worker_energy":
        worker = path.parent / "component_0_step_0_plus/result.json"
        result = json.loads(worker.read_text())
        result["energy_kcal_per_mol"] += 1
        _write(worker, result)
    else:
        _write(
            path.parent / "component_0_step_0_plus/process_status.json",
            {"exit_status": 2},
        )
    _write(path, check)
    with pytest.raises(ValueError):
        _prepare(prepared_inputs)


@pytest.mark.parametrize("run_requested", [False, True])
def test_cli_runs_engine_only_with_explicit_flag(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, run_requested: bool
) -> None:
    """Preparing through the public CLI must never implicitly execute Sander."""
    from scripts import run_nyl12_es_dftb_relaxation as cli

    executed = []
    monkeypatch.setattr(cli, "prepare_es_relaxation", lambda *args: {"prepared": True})
    monkeypatch.setattr(cli, "run_prepared_relaxation", executed.append)
    arguments = ["run_nyl12_es_dftb_relaxation.py"]
    for option in ["engine-config", "force-check", "region-config", "sander", "output"]:
        arguments.extend(["--" + option, str(tmp_path / option)])
    if run_requested:
        arguments.append("--run")
    monkeypatch.setattr(cli.sys, "argv", arguments)
    cli.main()
    assert executed == ([{"prepared": True}] if run_requested else [])
