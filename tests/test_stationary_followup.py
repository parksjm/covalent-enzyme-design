"""Exercise diagnostic control boundaries without launching a native engine."""

from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import stationary_candidate_branch as branch
from nyl12_covalent import stationary_followup as followup
from nyl12_covalent.stationary_native import digest


class HarmonicState:
    """Small analytic model embedded in the declared active array shape."""

    def __init__(self, source: Path) -> None:
        self.xyz = np.zeros((4165, 3))
        self.active = np.ones(4165, dtype=bool)
        self.indices = np.arange(4165)
        self.inv_sqrt = np.ones(12495)
        self.mdin = (
            "maxcyc=1, qmcharge=-2, dftb_telec=0.0, "
            "optalg='LBFGS', dimer=.false., maxcycle=1, tol=1.d-4, tole=1.d-5"
        )
        self.restart = object()
        self.topology = source / "topology"
        self.sander = Path("/fake/sander")

    def displaced(self, column: int, delta: float) -> np.ndarray:
        result = self.xyz.copy()
        result.reshape(-1)[column] += delta
        return result

    def gradient(self, xyz: np.ndarray, evidence: Path) -> tuple[float, np.ndarray]:
        diagonal = np.ones(12495)
        diagonal[:3] = [-0.23, 3.4e-5, 3.8e-5]
        gradient = xyz * diagonal.reshape(-1, 3)
        evidence.mkdir(parents=True)
        np.save(evidence / "gradient.npy", gradient)
        return float(np.sum(xyz * gradient) / 2), gradient


@pytest.fixture
def harmonic(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Provide immutable hashed columns, analytic state and in-process futures."""
    source, output = tmp_path / "source", tmp_path / "output"
    hessian = source / "ts_hessian"
    hessian.mkdir(parents=True)
    output.mkdir()
    modes = np.zeros((12495, 3))
    modes[:3, :3] = np.eye(3)
    np.save(hessian / "lowest_eigenvectors.npy", modes)
    np.save(hessian / "eigenvalues.npy", [-0.23, 3.4e-5, 3.8e-5])
    for col in (6586, 6587):
        folder = hessian / "columns" / f"column_{col:05d}"
        folder.mkdir(parents=True)
        vector = np.zeros(12495)
        vector[col] = 1.0
        np.save(folder / "column.npy", vector)
        (folder / "column.json").write_text(
            json.dumps({"step": 0.00025, "sha256": digest(folder / "column.npy")})
        )
    monkeypatch.setattr(followup, "NativeState", HarmonicState)
    monkeypatch.setattr(branch, "NativeState", HarmonicState)
    monkeypatch.setattr(followup, "ProcessPoolExecutor", ThreadPoolExecutor)
    monkeypatch.setattr(followup.os, "nice", lambda value: None)
    return source, output


def test_curvatures_and_independent_branch_failure(
    harmonic: tuple[Path, Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A failed branch preserves the other result and never asserts TS proof."""

    def fake_branch(args: tuple) -> dict:
        if args[2] == "minus":
            raise RuntimeError("one independent native failure")
        return {"geometry_class": "acyl_enzyme", "verified_ts": False}

    monkeypatch.setattr(followup, "candidate_branch", fake_branch)
    source, output = harmonic
    followup.run(source, output)
    curvature = json.loads((output / "curvature_diagnostic.json").read_text())
    assert curvature["candidate_descents_ready"]
    assert not curvature["verified_ts"]
    assert len(curvature["curvatures"]) == 9
    assert curvature["curvatures"][0]["gradient_curvature"] == pytest.approx(-0.23)
    branches = json.loads(
        (output / "candidate_connectivity_diagnostic.json").read_text()
    )
    assert "error" in branches["branches"]["minus"]
    assert branches["branches"]["plus"]["geometry_class"] == "acyl_enzyme"


def test_probe_failure_preserves_successful_evidence(
    harmonic: tuple[Path, Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Failed native evidence blocks descents while independent probes complete."""
    original = followup.probe

    def failed_probe(args: tuple) -> dict:
        if args[2:] == ("mode", 0, 0.01, -1):
            raise RuntimeError("SCC rejected")
        return original(args)

    monkeypatch.setattr(followup, "probe", failed_probe)
    source, output = harmonic
    with pytest.raises(RuntimeError, match="Native probes failed"):
        followup.run(source, output)
    receipt = json.loads((output / "probe_receipts.json").read_text())
    assert len(receipt["records"]) == 25
    assert len(receipt["failures"]) == 1
    assert not (output / "candidate_descents").exists()


def test_corrupt_original_column_rejected(harmonic: tuple[Path, Path]) -> None:
    """Refined probes cannot silently reuse a column with a changed hash."""
    source, output = harmonic
    receipt = source / "ts_hessian/columns/column_06586/column.json"
    receipt.write_text(json.dumps({"step": 0.00025, "sha256": "changed"}))
    with pytest.raises(ValueError, match="hash/step"):
        followup.run(source, output)


@pytest.mark.parametrize("outcome", ["success", "electronic", "process", "exception"])
def test_candidate_branch_native_boundaries(
    harmonic: tuple[Path, Path],
    monkeypatch: pytest.MonkeyPatch,
    outcome: str,
) -> None:
    """Record process, electronic, force and geometric outcomes separately."""
    source, output = harmonic
    monkeypatch.setattr(
        branch, "write_restart", lambda path, xyz, rst: path.write_text("restart")
    )
    final = np.zeros((19929, 3))
    final[13915] = [0, 0, 0]
    final[19915] = [1.4, 0, 0]
    final[19928] = [4, 0, 0]
    final[13916] = [3, 0, 0]
    monkeypatch.setattr(
        branch.Rst7, "open", lambda path: SimpleNamespace(coordinates=final)
    )
    monkeypatch.setattr(
        HarmonicState,
        "gradient",
        lambda self, xyz, evidence: (-1.0, np.zeros((4165, 3))),
    )

    def native(command: list, **kwargs: object) -> SimpleNamespace:
        if outcome == "exception":
            raise OSError("native unavailable")
        directory = kwargs["cwd"]
        mdin = (directory / "minimize.mdin").read_text()
        assert "maxcycle=20000" in mdin
        assert "dimer=.false." in mdin
        assert "optalg='LBFGS'" in mdin
        output = "SCC not converged"
        if outcome != "electronic":
            output = "\n".join(
                [
                    "Optimisation converged",
                    "QMMM SCC-DFTB:      Telec :     0.000K",
                    " QMMM SCC-DFTB: SCC-DFTB FOR STEP 1 DID NOT CONVERGE AFTER 70 cycles.",
                    " QMMM SCC-DFTB: SCC-DFTB for step 1 converged in 14 cycles.",
                ]
            )
        (directory / "minimize.mdout").write_text(output)
        (directory / "optimized.ncrst").touch()
        return SimpleNamespace(returncode=1 if outcome == "process" else 0)

    monkeypatch.setattr(branch.subprocess, "run", native)
    if outcome == "exception":
        with pytest.raises(OSError):
            branch.candidate_branch((source, output, "minus", 0.0))
        assert (output / "candidate_descents/minus/failure.json").exists()
    else:
        result = branch.candidate_branch((source, output, "plus", 0.0))
        assert not result["verified_ts"]
        if outcome == "success":
            assert result["force_converged"]
            assert result["geometry_class"] == "acyl_enzyme"
            assert result["electronic_gate"]
            assert result["electronic_history"]["recovered_trial_count"] == 1
            assert result["electronic_history"]["verified_scc_history"]
        else:
            assert "force_converged" not in result
