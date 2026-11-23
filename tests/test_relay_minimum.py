"""Bounded controller and numerical minimum acceptance regressions."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import relay_minimum as gates
from nyl12_covalent import relay_minimum_run as run
from scripts import watch_relay_minimum as watch


def geometry(og: float = 1.6, cn: float = 1.5, co: float = 1.35) -> dict:
    return {
        "distances_A": {"C_OG": og, "C_NL": cn, "C_O": co},
        "proton_bonds": {"HG": ["NL"]},
        "carbon_angles_degrees": {"O_C_OG": 108.0},
    }


@pytest.mark.parametrize(
    "og,cn,co,decision",
    [
        (1.6, 1.5, 1.35, "retained"),
        (1.8, 1.5, 1.35, "ambiguous_boundary"),
        (1.6, 1.7, 1.35, "ambiguous_boundary"),
        (1.6, 1.5, 1.28, "ambiguous_boundary"),
        (2.5, 1.4, 1.25, "outside_screen"),
        (1.6, 1.5, 1.2805, "ambiguous_boundary"),
    ],
)
def test_geometry_no_proton_label_restriction(og, cn, co, decision):
    result = gates.geometry_screen(geometry(og, cn, co))
    assert result["decision"] == decision
    assert result["geometry"]["proton_bonds"]["HG"] == ["NL"]


@pytest.mark.parametrize("value", [-1, np.nan])
def test_geometry_invalid_margin(value):
    with pytest.raises(ValueError):
        gates.geometry_screen(geometry(), value)


def test_geometry_rejects_nonfinite():
    with pytest.raises(ValueError):
        gates.geometry_screen(geometry(np.nan))


def test_repeat_force_and_serialization_gates():
    g = np.full((3, 3), 1e-6)
    assert run.gradient_metrics([1.0, 1.0, 1.0], [g, g, g])["passed"]
    assert not run.gradient_metrics([1.0, 1.0, 1.0 + 1e-12], [g, g, g])["passed"]
    assert not run.gradient_metrics([1.0 + 1e-7, 1.0, 1.0], [g, g, g])["passed"]
    assert not run.gradient_metrics([1.0, 1.0, 1.0], [g + 1e-5, g, g])["passed"]
    assert not run.gradient_metrics([1.0, 1.0, 1.0], [g, g, g + 1e-14])["passed"]
    assert not run.gradient_metrics([1.0, 1.0, 1.0], [g * 100, g * 100, g * 100])[
        "passed"
    ]
    for es, gs in [
        ([], []),
        ([1.0, 1.0, np.nan], [g, g, g]),
        ([1.0, 1.0, 1.0], [g, g, np.ones(2)]),
    ]:
        with pytest.raises(ValueError):
            run.gradient_metrics(es, gs)


def test_column_stability_explicit_engineering_gate():
    column = np.ones(10)
    assert run.column_metrics(column, column * (1 + 1e-5))["passed"]
    assert not run.column_metrics(column, column * 1.01)["passed"]
    for a, b in [
        (column, np.zeros(10)),
        (column, np.ones(2)),
        (column, column * np.nan),
    ]:
        with pytest.raises(ValueError):
            run.column_metrics(a, b)


def test_freeze_hashes_copy_and_reject_mutation(tmp_path):
    seed = tmp_path / "seed"
    reference = tmp_path / "ref" / "validated"
    output = tmp_path / "out"
    for p in (seed / "inputs", reference / "library_build", output):
        p.mkdir(parents=True)
    for n in ("system.prmtop", "dimer.mdin", "scientific_acceptance.json"):
        (seed / "inputs" / n).write_text(n)
    for n in ("optimized.ncrst", "scientific_acceptance.json", "source_identity.json"):
        (seed / n).write_text(n)
    (reference / "library_build/libdlfind.so").write_bytes(b"exact_library")
    (reference.parent / "build_receipt.json").write_text("{}")
    frozen = gates.freeze(seed, reference, output)
    gates.verify_frozen(frozen)
    assert (frozen / "inputs/optimized.ncrst").read_bytes() == (
        seed / "optimized.ncrst"
    ).read_bytes()
    (frozen / "inputs/optimized.ncrst").write_text("changed")
    with pytest.raises(ValueError, match="Frozen input changed"):
        gates.verify_frozen(frozen)


def test_job_records_do_not_infer_science_from_process(tmp_path):
    for name in gates.NAMES:
        seed = tmp_path / name
        (seed / "job").mkdir(parents=True)
        (seed / "job/status.json").write_text('{"exit_status":0}')
    assert all(r["scientific"] is None for r in gates.job_records(tmp_path).values())


def test_missing_amberhome_fails_explicitly(monkeypatch, tmp_path):
    monkeypatch.delenv("AMBERHOME", raising=False)
    with pytest.raises(KeyError):
        gates.prerequisites(tmp_path, tmp_path)
    monkeypatch.setenv("AMBERHOME", "/wrong")
    with pytest.raises(ValueError):
        gates.prerequisites(tmp_path, tmp_path)


def test_journal_requires_success(monkeypatch, tmp_path):
    monkeypatch.setenv("NYL12_REPOSITORY", str(tmp_path))
    calls = []
    monkeypatch.setattr(
        gates.subprocess, "run", lambda command, check: calls.append((command, check))
    )
    gates.journal("intent", "bounded stage", (tmp_path / "file",))
    assert calls[0][1] and "--file" in calls[0][0]


@pytest.mark.parametrize("passed", [True, False])
def test_full_hessian_expect_zero_and_keep_literal_result(
    monkeypatch, tmp_path, passed
):
    monkeypatch.setattr(run, "precision_preflight", lambda *args: {"passed": True})
    monkeypatch.setattr(run, "verify_frozen", lambda *args: None)
    monkeypatch.setattr(run, "verify_prerequisites", lambda *args: None)
    monkeypatch.setattr(run, "journal", lambda *args: None)
    calls = []

    def hessian(*args, **kwargs):
        calls.append((args, kwargs))
        return {"passed": passed, "negative": 0 if passed else 1, "unresolved": 0}

    monkeypatch.setattr(run, "full_hessian", hessian)
    result = run.verify_minimum(tmp_path / "frozen", tmp_path, 8)
    assert calls[0][1]["expected_negative"] == 0
    assert calls[0][0][3:] == (8, 0.0005)
    assert result["numerically_verified"] == passed
    assert "not a rigorous" in result["caveat"]


def test_failed_preflight_blocks_hessian(monkeypatch, tmp_path):
    monkeypatch.setattr(run, "precision_preflight", lambda *args: {"passed": False})
    assert not run.verify_minimum(tmp_path, tmp_path, 1)["numerically_verified"]
    with pytest.raises(ValueError):
        run.verify_minimum(tmp_path, tmp_path, 9)


def setup_watch(monkeypatch, tmp_path, records):
    monkeypatch.setattr(watch, "journal", lambda *args: None)
    monkeypatch.setattr(watch, "prerequisites", lambda *args: {"original": "sha"})
    monkeypatch.setattr(watch, "job_records", lambda *args: records)
    monkeypatch.setattr(
        watch, "NativeState", lambda *args: SimpleNamespace(xyz=np.zeros((1, 3)))
    )
    monkeypatch.setattr(watch, "geometry", lambda xyz: geometry())
    monkeypatch.setattr(watch, "proton_inventory", lambda *args: {})
    monkeypatch.setattr(watch, "resource_gate", lambda *args: {"passed": True})
    monkeypatch.setattr(watch, "freeze", lambda *args: tmp_path)
    return tmp_path / "campaign", tmp_path / "reference", tmp_path


def records(exit_status, scientific):
    return {
        n: {
            "process": {"exit_status": exit_status if n == "ti1_window033" else None},
            "scientific": scientific if n == "ti1_window033" else None,
        }
        for n in gates.NAMES
    }


def test_w033_triggers_while_other_jobs_run(monkeypatch, tmp_path):
    args = setup_watch(
        monkeypatch, tmp_path, records(0, {"accepted": True, "force_converged": True})
    )
    calls = []
    monkeypatch.setattr(watch, "run_owned", lambda *args: calls.append(args) or 0)
    watch.watch(*args, 120, 168, 8)
    assert len(calls) == 1 and "--verify" in calls[0][0]
    assert (
        json.loads((tmp_path / "process_acceptance.json").read_text())["exit_status"]
        == 0
    )


@pytest.mark.parametrize(
    "code,science",
    [(1, {"accepted": False}), (0, {"accepted": True, "force_converged": False})],
)
def test_no_restart_after_terminal_rejection(monkeypatch, tmp_path, code, science):
    args = setup_watch(monkeypatch, tmp_path, records(code, science))
    watch.watch(*args, 120, 168, 8)
    assert (
        json.loads((tmp_path / "scientific_acceptance.json").read_text())["status"]
        == "requires_review_or_continuation"
    )


def test_terminal_missing_record_is_error(monkeypatch, tmp_path):
    args = setup_watch(monkeypatch, tmp_path, records(0, None))
    with pytest.raises(RuntimeError, match="missing scientific"):
        watch.watch(*args, 120, 168, 8)


def test_ambiguous_geometry_blocks_derivatives(monkeypatch, tmp_path):
    args = setup_watch(
        monkeypatch, tmp_path, records(0, {"accepted": True, "force_converged": True})
    )
    monkeypatch.setattr(watch, "geometry", lambda xyz: geometry(1.8))
    watch.watch(*args, 120, 168, 8)
    assert (
        json.loads((tmp_path / "scientific_acceptance.json").read_text())["status"]
        == "ambiguous_boundary"
    )


def test_timeout_never_signals_original_jobs(monkeypatch, tmp_path):
    args = setup_watch(monkeypatch, tmp_path, records(None, None))
    ticks = iter([0, 2])
    monkeypatch.setattr(watch.time, "monotonic", lambda: next(ticks))
    with pytest.raises(TimeoutError):
        watch.watch(*args, 120, 0.0001, 8)


def test_nonzero_verifier_is_explicit_failure(monkeypatch, tmp_path):
    args = setup_watch(
        monkeypatch, tmp_path, records(0, {"accepted": True, "force_converged": True})
    )
    monkeypatch.setattr(watch, "run_owned", lambda *args: 124)
    with pytest.raises(RuntimeError, match="124"):
        watch.watch(*args, 120, 168, 8)


def test_live_records_ignore_partial_scientific_json(tmp_path):
    for name in gates.NAMES:
        seed = tmp_path / name
        (seed / "job").mkdir(parents=True)
        (seed / "job/status.json").write_text('{"exit_status":null}')
        (seed / "scientific_acceptance.json").write_text("{partial")
    assert all(r["scientific"] is None for r in gates.job_records(tmp_path).values())


def test_all_original_prerequisites_pinned(monkeypatch, tmp_path):
    campaign = tmp_path / "campaign"
    ref = tmp_path / "validated"
    amber = tmp_path / "amber"
    seed = campaign / "ti1_window033"
    for directory in (
        seed / "inputs",
        ref / "inputs",
        ref / "library_build",
        amber / "bin",
        amber / "dat/slko/3ob-3-1",
        campaign / "parameters_3ob",
    ):
        directory.mkdir(parents=True)
    monkeypatch.setattr(gates, "AMBER", amber)
    monkeypatch.setenv("AMBERHOME", str(amber))
    exe = amber / "bin/sander"
    exe.write_text("native")
    exe.chmod(0o700)
    lib = ref / "library_build/libdlfind.so"
    lib.write_text("precise")
    build = {
        "change": "write_xyz format only: f15.7 to es25.16",
        "library_sha256": gates.digest(lib),
    }
    (ref.parent / "build_receipt.json").write_text(json.dumps(build))
    params = {}
    from itertools import product

    for a, b in product("HCNO", repeat=2):
        name = f"{a}-{b}.skf"
        for directory in (amber / "dat/slko/3ob-3-1", campaign / "parameters_3ob"):
            (directory / name).write_text(name)
        params[name] = gates.digest(campaign / "parameters_3ob" / name)
    prepared = campaign / "prepared.json"
    prepared.write_text(json.dumps({"parameter_sha256": params}))
    identity = {
        "active_atom_count": 4165,
        "qm_atom_count": 98,
        "qm_mm_link_count": 8,
        "qm_charge": -2,
    }
    (seed / "source_identity.json").write_text(json.dumps(identity))
    for name in ("system.prmtop", "dimer.mdin", "scientific_acceptance.json"):
        for directory in (seed / "inputs", ref / "inputs"):
            (directory / name).write_text(name)
    hashes = gates.prerequisites(campaign, ref)
    assert len(hashes) == 40
    (amber / "dat/slko/3ob-3-1/H-H.skf").write_text("changed")
    with pytest.raises(ValueError, match="parameter changed"):
        gates.prerequisites(campaign, ref)
    (amber / "dat/slko/3ob-3-1/H-H.skf").write_text("H-H.skf")
    (seed / "inputs/dimer.mdin").write_text("changed")
    with pytest.raises(ValueError, match="reference input"):
        gates.prerequisites(campaign, ref)
    (seed / "inputs/dimer.mdin").write_text("dimer.mdin")
    lib.write_text("unqualified")
    with pytest.raises(ValueError, match="library differs"):
        gates.prerequisites(campaign, ref)
    lib.write_text("precise")
    build["change"] = "physics changed"
    (ref.parent / "build_receipt.json").write_text(json.dumps(build))
    with pytest.raises(ValueError, match="Unqualified"):
        gates.prerequisites(campaign, ref)


def test_probe_columns_include_exact_reactive_and_hb_atoms(monkeypatch):
    atoms = [SimpleNamespace(name="C") for _ in range(19916)]
    atoms[13924].name = "HB"
    monkeypatch.setattr(gates, "AmberParm", lambda path: SimpleNamespace(atoms=atoms))
    monkeypatch.setattr(
        gates,
        "AmberMask",
        lambda *args: SimpleNamespace(Selection=lambda: [False] * 19916),
    )
    state = SimpleNamespace(
        topology="top",
        mdin_source="qmmask='@1'",
        active=np.ones(19916, dtype=bool),
        indices=np.arange(19916),
    )
    columns = gates.probe_columns(state)
    assert set(range(13924 * 3, 13924 * 3 + 3)) <= set(columns)
    assert set(range(19915 * 3, 19915 * 3 + 3)) <= set(columns)
    atoms[13924].name = "HA"
    with pytest.raises(ValueError, match="mapping changed"):
        gates.probe_columns(state)


def test_proton_inventory_preserves_transfer_alternatives(monkeypatch):
    h = SimpleNamespace(idx=0, name="HG", atomic_number=1, bond_partners=[])
    atoms = [
        h,
        SimpleNamespace(idx=1, atomic_number=7),
        SimpleNamespace(idx=2, atomic_number=8),
    ]
    monkeypatch.setattr(gates, "AmberParm", lambda path: SimpleNamespace(atoms=atoms))
    monkeypatch.setattr(
        gates, "AmberMask", lambda *args: SimpleNamespace(Selection=lambda: [True] * 3)
    )
    state = SimpleNamespace(
        topology="top",
        mdin_source="qmmask='@1-3'",
        xyz=np.array([[0, 0, 0], [1, 0, 0], [1.21, 0, 0]]),
    )
    result = gates.proton_inventory(state)
    assert result["qm_hydrogen_count"] == 1
    assert result["protons"][0]["partners_within_1p2A"] == [2]
    assert result["protons"][0]["ambiguous_cutoff"]


@pytest.mark.parametrize("passing", [True, False])
def test_native_preflight_blocks_or_executes_columns(monkeypatch, tmp_path, passing):
    monkeypatch.setattr(run, "verify_frozen", lambda *args: None)
    monkeypatch.setattr(run, "verify_prerequisites", lambda *args: None)
    monkeypatch.setattr(run, "journal", lambda *args: None)
    monkeypatch.setattr(
        run, "proton_inventory", lambda *args: {"qm_hydrogen_count": 42}
    )
    monkeypatch.setattr(run, "probe_columns", lambda *args: [0])

    class State:
        def __init__(self, *args):
            self.indices = np.arange(4165)
            self.xyz = np.zeros((4165, 3))
            self.restart_path = tmp_path / "restart"

        def gradient(self, *args, **kwargs):
            return 1.0, np.full((4165, 3), 1e-6 if passing else 1e-3)

    monkeypatch.setattr(run, "NativeState", State)

    def worker(task):
        root, stage, restart, col, step = task
        folder = stage / "columns" / f"column_{col:05d}"
        folder.mkdir(parents=True)
        np.save(folder / "column.npy", np.ones(12495))
        return str(folder / "column.npy")

    class Executor:
        def __init__(self, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def map(self, func, tasks):
            return map(worker, tasks)

    monkeypatch.setattr(run, "ProcessPoolExecutor", Executor)
    result = run.precision_preflight(tmp_path, tmp_path, 8)
    assert result["passed"] == passing
    assert (tmp_path / "03_step_probes").exists() == passing


def test_original_prerequisites_cannot_change_between_stages(tmp_path):
    source = tmp_path / "original"
    source.write_text("original")
    (tmp_path / "prerequisites.json").write_text(
        json.dumps({str(source): gates.digest(source)})
    )
    run.verify_prerequisites(tmp_path)
    source.write_text("changed")
    with pytest.raises(ValueError):
        run.verify_prerequisites(tmp_path)


def test_owned_verifier_normal_process(monkeypatch, tmp_path):
    class Process:
        pid = 1234

        def wait(self, timeout=None):
            return 0

    monkeypatch.setattr(watch.subprocess, "Popen", lambda *args, **kwargs: Process())
    assert watch.run_owned(["verification"], tmp_path, 3) == 0


def test_resource_headroom_actual_readonly(tmp_path):
    result = watch.resource_gate(tmp_path, 8)
    assert result["workers"] == 8 and "available_memory_bytes" in result


@pytest.mark.parametrize(
    "hours,interval", [(np.inf, 120), (168, np.nan), (0, 120), (168, 0)]
)
def test_controller_rejects_unbounded_or_invalid_timing(tmp_path, hours, interval):
    with pytest.raises(ValueError):
        watch.watch(tmp_path, tmp_path, tmp_path, interval, hours, 8)
