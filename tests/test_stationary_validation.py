"""Scientific correctness and rejection tests for complete Hessian validation."""

from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import stationary_endpoints as end
from nyl12_covalent import stationary_hessian as hs
from nyl12_covalent import stationary_native as native
from nyl12_covalent import stationary_precision as precision


def test_cartesian_difference_applies_both_mass_factors(tmp_path, monkeypatch):
    """A coupled analytic quadratic gives the exact mass-weighted column."""
    matrix = np.array([[4.0, 2.0, 0.0], [2.0, 8.0, 1.0], [0.0, 1.0, 3.0]])
    state = SimpleNamespace(inv_sqrt=np.array([0.5, 0.5, 0.5]))
    state.displaced = lambda column, delta: np.eye(3)[column] * delta
    state.gradient = lambda xyz, evidence: (0.0, matrix @ xyz)
    actual = precision.column_difference(state, 1, 0.002, tmp_path)
    np.testing.assert_allclose(actual, matrix[:, 1] / 4)
    monkeypatch.setattr(native, "cached_state", lambda root, restart: state)
    restart = tmp_path / "restart"
    restart.write_text("immutable")
    stage = tmp_path / "hessian"
    output = Path(native.compute_column(tmp_path, stage, restart, 1, 0.002))
    np.testing.assert_allclose(np.load(output), actual)
    assert native.compute_column(tmp_path, stage, restart, 1, 0.002) == str(output)
    with pytest.raises(ValueError, match="digest/step"):
        native.compute_column(tmp_path, stage, restart, 1, 0.001)


@pytest.mark.parametrize(
    "values,expected,passed,unresolved",
    [
        ([-2.0, 0.1, 1.0], 1, True, 0),
        ([-2.0, -0.1, 1.0], 1, False, 0),
        ([-2.0, 1e-8, 1.0], 1, False, 1),
        ([0.1, 0.2, 1.0], 0, True, 0),
    ],
)
def test_inertia_requires_all_signs_resolved(values, expected, passed, unresolved):
    record = hs.validate_inertia(np.array(values), np.array([1e-12]), expected, 1e-6)
    assert record["passed"] is passed
    assert record["unresolved"] == unresolved


def test_inertia_rejects_bad_residual_and_nonfinite():
    assert not hs.validate_inertia(np.array([-1.0, 1.0]), np.array([0.01]), 1, 0.001)[
        "passed"
    ]
    with pytest.raises(ValueError):
        hs.validate_inertia(np.array([np.nan]), np.array([0.0]), 0, 0.001)


def test_complete_hessian_recovers_known_inertia(tmp_path, monkeypatch):
    root = tmp_path / "root"
    (root / "library_build").mkdir(parents=True)
    (root / "library_build/libdlfind.so").write_text("writer")
    restart = root / "restart"
    restart.write_text("coordinate")
    stage = root / "hessian"
    rng = np.random.default_rng(42)
    rotation, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    matrix = rotation @ np.diag([-2.0, 0.1, 0.2, 0.5, 1.0, 3.0]) @ rotation.T
    monkeypatch.setattr(
        hs, "NativeState", lambda *args: SimpleNamespace(indices=[0, 1])
    )
    monkeypatch.setattr(hs, "ProcessPoolExecutor", ThreadPoolExecutor)

    def column(args):
        index = args[3]
        out = stage / "columns" / f"column_{index:05d}"
        out.mkdir(parents=True)
        np.save(out / "column.npy", matrix[:, index])
        return str(out / "column.npy")

    monkeypatch.setattr(hs, "worker_column", column)
    result = hs.full_hessian(root, restart, stage, 2, 0.002, 1)
    assert result["passed"] is True
    assert result["dimension"] == 6
    np.testing.assert_allclose(
        np.load(stage / "eigenvalues.npy"), [-2.0, 0.1, 0.2, 0.5, 1.0, 3.0]
    )
    with pytest.raises(ValueError, match="specification"):
        hs.full_hessian(root, restart, stage, 2, 0.001, 1)


def test_precision_gate_rejects_old_new_mismatch(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    root.mkdir()
    (root / "inputs").mkdir()
    (root / "inputs/native_candidate.mdout").write_text("Final converged energy:")

    class FakeState:
        active = np.ones(4165, dtype=bool)
        indices = np.arange(4165)
        xyz = np.zeros((4165, 3))

        def gradient(self, xyz, evidence, precision=True):
            return (-1.0, np.full((4165, 3), 1e-6 if precision else 0.0))

    monkeypatch.setattr(precision, "NativeState", lambda root: FakeState())
    monkeypatch.setattr(
        precision, "read_xyz_frames", lambda *args: [np.ones((4165, 3))]
    )
    monkeypatch.setattr(precision, "column_difference", lambda *args: np.ones(3))
    assert precision.precision_gate(root, 0.002)["passed"] is False


@pytest.mark.parametrize(
    "distances,label",
    [
        ({"C_N": 1.35, "OG_C": 3.0, "OG_H": 1.0, "N_H": 2.0}, "reactant_ES"),
        ({"C_N": 3.0, "OG_C": 1.35, "OG_H": 2.0, "N_H": 1.0}, "acyl_enzyme"),
        ({"C_N": 1.35, "OG_C": 1.35, "OG_H": 2.0, "N_H": 2.0}, "unassigned"),
    ],
)
def test_chemical_labels_require_proton_and_bond_inventory(distances, label):
    assert end.classify_distances(distances) == label


@pytest.mark.parametrize("frame_shift", [0.0, 0.01])
def test_native_gradient_archives_lossless_evidence(tmp_path, monkeypatch, frame_shift):
    state = native.NativeState.__new__(native.NativeState)
    state.root = tmp_path
    state.active = np.array([True, False])
    state.indices = np.array([0])
    state.topology = tmp_path / "topology"
    state.restart = None
    state.mdin = "native matched"
    state.sander = tmp_path / "sander"
    state.scale = 1.0

    class Evaluator:
        def __init__(self, topology, restart, active, mdin, sander, branch, scale):
            self.branch = branch
            self.mdin = mdin
            self.evidence = [{"gradient": "exact"}]

        def __call__(self, xyz):
            for name, text in {
                "eval_00001.mdout": "DFTBESCF= -3\nRun   done\nWorking directory: /tmp/nyl12_gradient_9fnanpmc",
                "eval_00001.path_force.xyz": "1\n\nH 1.0 2.0 3.0\n",
                "eval_00001.mdin": self.mdin,
                "eval_00001.path.xyz": f"1\n\nH {frame_shift} 0 0\n",
                "eval_00001.log": "",
            }.items():
                (self.branch / name).write_text(text)
            return -1.0, np.ones((2, 3))

    monkeypatch.setattr(native, "NativeSanderEvaluator", Evaluator)
    if frame_shift:
        with pytest.raises(RuntimeError, match="frame identity"):
            state.gradient(np.zeros((2, 3)), tmp_path / "evidence")
        assert (tmp_path / "evidence/eval_00001.path.xyz.gz").exists()
        return
    energy, gradient = state.gradient(np.zeros((2, 3)), tmp_path / "evidence")
    assert energy == -1.0
    cached_energy, cached_gradient = state.gradient(
        np.zeros((2, 3)), tmp_path / "evidence"
    )
    assert cached_energy == energy
    np.testing.assert_array_equal(cached_gradient, gradient)
    np.testing.assert_allclose(gradient, np.ones((1, 3)) / 627.5094740631)
    receipt = json.loads((tmp_path / "evidence/receipt.json").read_text())
    assert receipt["serialization"] == "es25.16"
    assert receipt["electronic_history"]["passed"]
    assert not receipt["electronic_history"]["verified_scc_history"]
    with gzip.open(tmp_path / "evidence/eval_00001.path_force.xyz.gz", "rt") as stream:
        assert stream.read() == "1\n\nH 1.0 2.0 3.0\n"
    assert (
        receipt["gradient_sha256"]
        == hashlib.sha256((tmp_path / "evidence/gradient.npy").read_bytes()).hexdigest()
    )


def test_native_state_preserves_full_active_mask_and_coordinate_displacements(
    tmp_path, monkeypatch
):
    (tmp_path / "inputs").mkdir()
    (tmp_path / "inputs/dimer.mdin").write_text("active=':1', dftb_telec=0.0")
    monkeypatch.setattr(
        native.Rst7,
        "open",
        lambda path: SimpleNamespace(coordinates=np.zeros((4165, 3))),
    )
    monkeypatch.setattr(
        native,
        "AmberParm",
        lambda path: SimpleNamespace(atoms=[SimpleNamespace(mass=4.0)] * 4165),
    )
    monkeypatch.setattr(
        native,
        "AmberMask",
        lambda *args: SimpleNamespace(Selection=lambda: [True] * 4165),
    )
    monkeypatch.setattr(native, "load_calibration", lambda path: (1.0, {}))
    monkeypatch.setattr(native, "prepare_native_input", lambda source: source)
    state = native.NativeState(tmp_path)
    assert len(state.indices) == 4165
    np.testing.assert_allclose(state.inv_sqrt, 0.5)
    displaced = state.displaced(4, 0.002)
    assert displaced[1, 1] == 0.002
    assert np.count_nonzero(displaced) == 1
    with pytest.raises(ValueError):
        state.displaced(12495, 0.002)
    (tmp_path / "inputs/dimer.mdin").write_text("dftb_telec=0.0")
    with pytest.raises(ValueError, match="active"):
        native.NativeState(tmp_path)


def test_both_verified_mode_directions_reach_distinct_chemical_endpoints(
    tmp_path, monkeypatch
):
    root = tmp_path / "root"
    root.mkdir()
    xyz = np.zeros((20000, 3))
    state = SimpleNamespace(
        xyz=xyz,
        active=np.ones(20000, dtype=bool),
        inv_sqrt=np.ones(60000),
        restart=None,
        mdin="maxcyc=1 maxcycle=1, tol=1.d-4, tole=1.d-5",
        topology=root / "topology",
        sander=Path("sander"),
    )
    state.gradient = lambda xyz, evidence: (-2.0, np.zeros((20000, 3)))
    monkeypatch.setattr(end, "NativeState", lambda *args: state)
    mode = np.zeros((60000, 1))
    mode[0, 0] = 1.0
    monkeypatch.setattr(end.np, "load", lambda *args: mode)
    monkeypatch.setattr(
        end, "write_restart", lambda path, *args: path.write_text("seed")
    )

    def minimize(command, cwd, **kwargs):
        (cwd / "minimize.mdout").write_text("Optimisation converged")
        (cwd / "optimized.ncrst").write_text("endpoint")

    monkeypatch.setattr(end.subprocess, "run", minimize)

    def read_restart(path):
        coords = xyz.copy()
        if Path(path).parent.name == "minus":
            coords[13916, 0], coords[19915, 0], coords[19928, 0] = 1.0, 3.0, 4.35
        else:
            coords[13916, 0], coords[19915, 0], coords[19928, 0] = 3.0, 1.35, 4.0
        return SimpleNamespace(coordinates=coords)

    monkeypatch.setattr(end.Rst7, "open", read_restart)
    points = end.endpoints(root, -1.0)
    assert len(points) == 2
    assert json.loads((root / "endpoints/connectivity_acceptance.json").read_text())[
        "passed"
    ]


def test_worker_does_not_change_existing_low_priority(monkeypatch):
    monkeypatch.setattr(hs.os, "getpriority", lambda *args: 10)
    monkeypatch.setattr(hs, "compute_column", lambda *args: "column")
    assert (
        hs.worker_column((Path("root"), Path("stage"), Path("restart"), 0, 0.002))
        == "column"
    )


def test_native_scc_warning_is_rejected_and_preserved(tmp_path, monkeypatch):
    state = native.NativeState.__new__(native.NativeState)
    for name, value in {
        "root": tmp_path,
        "active": np.array([True]),
        "topology": tmp_path / "topology",
        "restart": None,
        "mdin": "input",
        "sander": tmp_path / "sander",
        "scale": 1.0,
    }.items():
        setattr(state, name, value)

    class FailedElectronic:
        def __init__(self, *args):
            self.branch = args[-2]
            self.evidence = []

        def __call__(self, xyz):
            (self.branch / "eval_00001.mdout").write_text(
                "DFTBESCF= -3\nRun   done\n QMMM SCC-DFTB: Convergence could not be achieved in this step."
            )
            return -1.0, np.zeros((1, 3))

    monkeypatch.setattr(native, "NativeSanderEvaluator", FailedElectronic)
    with pytest.raises(RuntimeError, match="electronic"):
        state.gradient(np.zeros((1, 3)), tmp_path / "failed")
    assert (tmp_path / "failed/eval_00001.mdout.gz").exists()
    assert not (tmp_path / "failed/receipt.json").exists()


def test_failed_native_process_preserves_logs_and_short_topology(tmp_path, monkeypatch):
    state = native.NativeState.__new__(native.NativeState)
    for name, value in {
        "root": tmp_path,
        "active": np.array([True]),
        "topology": tmp_path / "long_topology.prmtop",
        "restart": None,
        "mdin": "input",
        "sander": tmp_path / "sander",
        "scale": 1.0,
    }.items():
        setattr(state, name, value)

    class FailedProcess:
        def __init__(self, *args):
            self.branch = args[-2]
            assert args[0].name == "topology.prmtop"
            assert args[0].is_symlink()

        def __call__(self, xyz):
            (self.branch / "eval_00001.log").write_text("native process failure reason")
            raise RuntimeError("native exit1")

    monkeypatch.setattr(native, "NativeSanderEvaluator", FailedProcess)
    with pytest.raises(RuntimeError, match="native exit1"):
        state.gradient(np.zeros((1, 3)), tmp_path / "failed")
    assert (
        json.loads((tmp_path / "failed/failure.json").read_text())["accepted"] is False
    )
    with gzip.open(tmp_path / "failed/eval_00001.log.gz", "rt") as stream:
        assert stream.read() == "native process failure reason"
