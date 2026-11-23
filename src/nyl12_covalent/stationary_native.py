"""Native gradients with reconstructible inputs and compact lossless evidence."""

from __future__ import annotations

import gzip
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any

import numpy as np

from scripts.prepare_dimer_continuation import read_xyz_frames
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from nyl12_covalent.dftb_scc_validation import audit_dftb_scc
from nyl12_covalent.native_sander_cli import (
    HARTREE_TO_KCAL_PER_MOL,
    NativeSanderEvaluator,
    load_calibration,
    prepare_native_input,
)


def digest(path: Path) -> str:
    """Return the SHA-256 digest of an existing artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Atomically publish a finite JSON checkpoint."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def archive_native_outputs(scratch: Path, evidence: Path) -> None:
    """Losslessly retain native text evidence, including failed process output."""
    for pattern in ("*.mdout", "*.log", "*.path_force.xyz", "*.path.xyz", "*.mdin"):
        for source in scratch.glob(pattern):
            target = evidence / (source.name + ".gz")
            with gzip.open(target, "wb") as stream:
                stream.write(source.read_bytes())
            with gzip.open(target, "rb") as stream:
                if hashlib.sha256(stream.read()).hexdigest() != digest(source):
                    raise RuntimeError("Lossless evidence compression failed")


class NativeState:
    """Immutable QM/MM state with the exact declared 4165-atom active shell."""

    def __init__(self, root: Path, restart: Path | None = None) -> None:
        self.root = root
        self.topology = root / "inputs/system.prmtop"
        self.restart_path = restart or root / "inputs/optimized.ncrst"
        self.restart = Rst7.open(str(self.restart_path))
        self.xyz = np.asarray(self.restart.coordinates).reshape((-1, 3))
        self.mdin_source = (root / "inputs/dimer.mdin").read_text()
        matches = re.findall(r"\bactive\s*=\s*'([^']+)'", self.mdin_source, re.I)
        if len(matches) != 1:
            raise ValueError("Require exactly one immutable active atom mask")
        parm = AmberParm(str(self.topology))
        self.active = np.asarray(AmberMask(parm, matches[0]).Selection(), dtype=bool)
        if self.active.sum() != 4165 or len(parm.atoms) != len(self.xyz):
            raise ValueError("Active region or topology/restart atom count changed")
        self.indices = np.flatnonzero(self.active)
        self.masses = np.asarray([a.mass for a in parm.atoms])[self.active]
        self.inv_sqrt = np.repeat(1 / np.sqrt(self.masses), 3)
        self.scale, _ = load_calibration(root / "inputs/scientific_acceptance.json")
        self.mdin = prepare_native_input(self.mdin_source)
        if not re.search(r"dftb_telec\s*=\s*0\.0", self.mdin):
            raise ValueError("Expected zero electronic temperature")
        self.sander = Path("/data/jpc/code/ambertools26/venv/bin/sander")

    def gradient(
        self, xyz: np.ndarray, evidence: Path, *, precision: bool = True
    ) -> tuple[float, np.ndarray]:
        """Return analytic native gradients and archive exact compact evidence.

        Args:
            xyz: Full coordinates reconstructed from retained state/displacement.
            evidence: Unique evidence directory. Full restart copies are scratch.
            precision: Select the isolated serialization-only library.

        Returns:
            Energy in Hartree and active Cartesian gradient in Hartree/Angstrom.
        """
        evidence.mkdir(parents=True, exist_ok=True)
        if (evidence / "receipt.json").exists():
            cached = json.loads((evidence / "receipt.json").read_text())
            expected_serialization = "es25.16" if precision else "f15.7"
            with gzip.open(next(evidence.glob("*.mdin.gz")), "rt") as stream:
                if stream.read() != self.mdin:
                    raise ValueError("Existing native gradient uses a different MDIN")
            if (
                cached["coordinates_sha256"]
                != hashlib.sha256(xyz.tobytes()).hexdigest()
                or cached["gradient_sha256"] != digest(evidence / "gradient.npy")
                or cached["serialization"] != expected_serialization
            ):
                raise ValueError(
                    "Existing native gradient does not match coordinates/hash/precision"
                )
            return float(cached["energy_hartree"]), np.load(evidence / "gradient.npy")
        if any(evidence.iterdir()):
            raise FileExistsError(
                "Incomplete native evidence must remain preserved in its attempt"
            )
        old = os.environ.get("LD_LIBRARY_PATH")
        amberlib = "/data/jpc/code/ambertools26/venv/lib"
        os.environ["LD_LIBRARY_PATH"] = (
            f"{self.root / 'library_build'}:{amberlib}" if precision else amberlib
        )
        try:
            with tempfile.TemporaryDirectory(prefix="nyl12_gradient_") as scratch:
                short_topology = Path(scratch) / "topology.prmtop"
                short_topology.symlink_to(self.topology)
                evaluator = NativeSanderEvaluator(
                    short_topology,
                    self.restart,
                    self.active,
                    self.mdin,
                    self.sander,
                    Path(scratch),
                    self.scale,
                )
                try:
                    energy, full_gradient = evaluator(xyz)
                except Exception as error:
                    archive_native_outputs(Path(scratch), evidence)
                    write_json(
                        evidence / "failure.json",
                        {"error": str(error), "accepted": False},
                    )
                    raise
                mdout = next(Path(scratch).glob("*.mdout"))
                text = mdout.read_text()
                electronic_history = audit_dftb_scc(text)
                archive_native_outputs(Path(scratch), evidence)
                if (
                    not electronic_history["passed"]
                    or "DFTBESCF=" not in text
                    or "Run   done" not in text
                ):
                    raise RuntimeError(
                        "Native electronic/termination gate failed: "
                        f"{electronic_history}"
                    )
                coordinate_frames = read_xyz_frames(
                    next(Path(scratch).glob("*.path.xyz")), len(self.indices)
                )
                force_frames = read_xyz_frames(
                    next(Path(scratch).glob("*.path_force.xyz")), len(self.indices)
                )
                frame_error = float(
                    np.max(
                        np.abs(
                            coordinate_frames[-1] / 0.999999999068864 - xyz[self.active]
                        )
                    )
                )
                frame_tolerance = 1e-10 if precision else 5.01e-8
                if (
                    len(coordinate_frames) != 1
                    or len(force_frames) != 1
                    or frame_error > frame_tolerance
                ):
                    write_json(
                        evidence / "frame_failure.json",
                        {
                            "coordinate_frames": len(coordinate_frames),
                            "force_frames": len(force_frames),
                            "error_angstrom": frame_error,
                        },
                    )
                    raise RuntimeError(
                        "Native evaluated-coordinate/force frame identity gate failed"
                    )
                gradient = full_gradient[self.active] / HARTREE_TO_KCAL_PER_MOL
                np.save(evidence / "gradient.npy", gradient, allow_pickle=False)
                receipt = {
                    "evaluated_coordinate_frame_identity": {
                        "passed": True,
                        "coordinate_frames": len(coordinate_frames),
                        "force_frames": len(force_frames),
                        "max_error_angstrom": frame_error,
                        "source_conversion_scale": 0.999999999068864,
                    },
                    "energy_hartree": energy,
                    "coordinates_sha256": hashlib.sha256(xyz.tobytes()).hexdigest(),
                    "gradient_sha256": digest(evidence / "gradient.npy"),
                    "electronic_gate": "finite DFTBESCF; accepted SCC history; native Run done",
                    "electronic_history": electronic_history,
                    "serialization": "es25.16" if precision else "f15.7",
                    "evidence": evaluator.evidence,
                    "scratch_policy": "Transient full-system restart copies discarded only after lossless native force/mdout/input archival and gradient hashing; coordinates reconstruct from immutable state and column displacement.",
                }
                write_json(evidence / "receipt.json", receipt)
                return energy, gradient
        finally:
            if old is None:
                os.environ.pop("LD_LIBRARY_PATH", None)
            else:
                os.environ["LD_LIBRARY_PATH"] = old

    def displaced(self, column: int, delta: float) -> np.ndarray:
        """Return a full geometry displaced in one active Cartesian coordinate."""
        if not 0 <= column < len(self.indices) * 3 or delta == 0:
            raise ValueError("Invalid active column or displacement")
        xyz = self.xyz.copy()
        xyz[self.indices[column // 3], column % 3] += delta
        return xyz


def compute_column(
    root: Path, stage: Path, restart: Path, column: int, step: float
) -> str:
    """Compute a restartable full Hessian column and retain both gradients.

    Args:
        root: Immutable calculation inputs and precision library.
        stage: New or validated continuation Hessian directory.
        restart: Exact stationary candidate geometry.
        column: Active Cartesian degree of freedom in topology order.
        step: Cartesian finite-difference displacement in Angstrom.

    Returns:
        Path of the completed mass-weighted Hessian column.
    """
    directory = stage / "columns" / f"column_{column:05d}"
    completed = directory / "column.npy"
    if completed.exists():
        receipt = json.loads((directory / "column.json").read_text())
        if receipt["sha256"] != digest(completed) or receipt["step"] != step:
            raise ValueError("Completed Hessian column digest/step mismatch")
        return str(completed)
    state = cached_state(root, restart)
    directory.mkdir(parents=True, exist_ok=False)
    _, plus = state.gradient(state.displaced(column, step), directory / "plus")
    _, minus = state.gradient(state.displaced(column, -step), directory / "minus")
    hessian = ((plus - minus) / (2 * step)).reshape(-1)
    hessian *= state.inv_sqrt * state.inv_sqrt[column]
    np.save(completed, hessian, allow_pickle=False)
    write_json(
        directory / "column.json",
        {
            "column": column,
            "step": step,
            "sha256": digest(completed),
            "restart_sha256": digest(restart),
            "units": "Hartree/(Angstrom^2 amu)",
        },
    )
    return str(completed)


@lru_cache(maxsize=2)
def cached_state(root: Path, restart: Path) -> NativeState:
    """Reuse immutable topology parsing within one worker process."""
    return NativeState(root, restart)
