"""Mocked OpenMM tests that do not create a native CUDA context."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest

from nyl12_covalent import design_openmm_control as worker


class FakeUnit:
    """Represent OpenMM units without performing unit arithmetic."""

    def __init__(self, name: str) -> None:
        self.name = name

    def __rmul__(self, value: float) -> FakeUnit:
        return FakeUnit(f"{value}*{self.name}")

    def __truediv__(self, other: FakeUnit) -> FakeUnit:
        return FakeUnit(f"{self.name}/{other.name}")


class FakeQuantity:
    """Return a stable array or scalar when converted to an OpenMM unit."""

    def __init__(self, value: Any) -> None:
        self.value = value

    def value_in_unit(self, _unit: FakeUnit) -> Any:
        return self.value


def _install_fake_openmm(
    monkeypatch: pytest.MonkeyPatch, *, missing_box: bool = False
) -> dict[str, Any]:
    """Install a no-GPU OpenMM facade and return captured platform settings."""
    captured: dict[str, Any] = {}
    units = SimpleNamespace(
        kilojoule_per_mole=FakeUnit("kJ/mol"),
        nanometer=FakeUnit("nm"),
        femtoseconds=FakeUnit("fs"),
    )

    class FakeContext:
        minimized = False

        def setPeriodicBoxVectors(self, *vectors: Any) -> None:
            captured["box_vectors"] = vectors

        def setPositions(self, positions: Any) -> None:
            captured["positions"] = positions

        def getState(self, **_kwargs: Any) -> Any:
            energy = 5.0 if self.minimized else 10.0
            magnitude = 0.1 if self.minimized else 1.0
            forces = np.full((2, 3), magnitude)
            positions = np.full((2, 3), 0.25)
            return SimpleNamespace(
                getPotentialEnergy=lambda: FakeQuantity(energy),
                getForces=lambda **_: FakeQuantity(forces),
                getPositions=lambda **_: FakeQuantity(positions),
            )

    context = FakeContext()

    class FakePrmtop:
        topology = object()

        def __init__(self, _path: str) -> None:
            pass

        def createSystem(self, **kwargs: Any) -> object:
            captured["system_options"] = kwargs
            return object()

    class FakeInpcrd:
        positions = object()
        boxVectors = None if missing_box else (1, 2, 3)

        def __init__(self, _path: str) -> None:
            pass

    class FakeSimulation:
        def __init__(
            self,
            _topology: Any,
            _system: Any,
            _integrator: Any,
            _platform: Any,
            properties: dict[str, str],
        ) -> None:
            captured["properties"] = properties
            self.context = context

    app = ModuleType("openmm.app")
    app.PME = object()
    app.AmberPrmtopFile = FakePrmtop
    app.AmberInpcrdFile = FakeInpcrd
    app.Simulation = FakeSimulation
    app.PDBFile = SimpleNamespace(
        writeFile=lambda _topology, _positions, stream, keepIds: stream.write(
            f"PDB keepIds={keepIds}\n"
        )
    )
    openmm = ModuleType("openmm")
    openmm.app = app
    openmm.unit = units
    openmm.VerletIntegrator = lambda _step: object()
    openmm.Platform = SimpleNamespace(
        getPlatformByName=lambda name: captured.setdefault("platform", name)
    )
    openmm.LocalEnergyMinimizer = SimpleNamespace(
        minimize=lambda _context, tolerance, max_iterations: (
            captured.update(
                minimized=True, tolerance=tolerance, max_iterations=max_iterations
            )
            or setattr(context, "minimized", True)
        )
    )
    monkeypatch.setitem(sys.modules, "openmm", openmm)
    monkeypatch.setitem(sys.modules, "openmm.app", app)
    monkeypatch.setitem(sys.modules, "openmm.unit", units)
    return captured


def _input_files(tmp_path: Path) -> tuple[Path, Path]:
    topology = tmp_path / "complex.prmtop"
    coordinates = tmp_path / "complex.rst7"
    topology.write_text("topology\n")
    coordinates.write_text("coordinates\n")
    return topology, coordinates


def test_run_minimization_writes_coordinates_and_acceptance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The worker records the explicit no-dynamics CUDA minimization settings."""
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    captured = _install_fake_openmm(monkeypatch)
    topology, coordinates = _input_files(tmp_path)
    output = tmp_path / "control"

    receipt = worker.run_minimization(topology, coordinates, output)

    assert np.load(output / "positions_nm.npy").shape == (2, 3)
    assert (output / "final.pdb").is_file()
    assert json.loads((output / "acceptance.json").read_text()) == receipt
    assert receipt["energy_decreased"] is True
    assert receipt["force_converged"] is True
    assert receipt["dynamics_performed"] is False
    assert receipt["final_rms_gradient_kcal_mol_angstrom"] == pytest.approx(0.1 / 41.84)
    assert captured["platform"] == "CUDA"
    assert captured["properties"] == {"DeviceIndex": "0", "Precision": "mixed"}
    assert captured["max_iterations"] == 10_000


def test_worker_requires_visible_device_zero(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A control cannot silently move to a different GPU."""
    topology, coordinates = _input_files(tmp_path)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "1")
    with pytest.raises(ValueError, match="CUDA_VISIBLE_DEVICES"):
        worker.run_minimization(topology, coordinates, tmp_path / "control")


def test_worker_rejects_restart_without_periodic_box(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The PME run requires box vectors from the matching restart."""
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    _install_fake_openmm(monkeypatch, missing_box=True)
    topology, coordinates = _input_files(tmp_path)
    with pytest.raises(ValueError, match="periodic box"):
        worker.run_minimization(topology, coordinates, tmp_path / "control")
