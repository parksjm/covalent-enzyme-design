"""Tests for the intrinsic proton-transfer barrier benchmark."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType

import numpy as np
import pytest


SOURCE = (
    Path(__file__).parents[1]
    / "qmmm/benchmarks/02_small_model_intrinsic_barriers/source"
)


def _load(path: Path, name: str) -> ModuleType:
    """Load one script as a module for focused tests."""
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _fixed_generator() -> ModuleType:
    """Load the preceding matched fixed-geometry generator."""
    path = (
        Path(__file__).parents[1]
        / "qmmm/benchmarks/01_small_model_proton_transfer/source/generate_benchmark.py"
    )
    return _load(path, "barrier_test_fixed_generator")


def _write_fchk(path: Path, numbers: list[int], coordinates_a: np.ndarray) -> None:
    """Write the small FCHK subset consumed by the converter."""
    coordinates_bohr = coordinates_a.reshape(-1) / 0.529177210903
    lines = [f"Atomic numbers                        I   N= {len(numbers):12d}"]
    lines.extend(
        " ".join(str(value) for value in numbers[index : index + 6])
        for index in range(0, len(numbers), 6)
    )
    lines.append(
        f"Current cartesian coordinates         R   N= {len(coordinates_bohr):12d}"
    )
    lines.extend(
        " ".join(f"{value: .8E}" for value in coordinates_bohr[index : index + 5])
        for index in range(0, len(coordinates_bohr), 5)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def prepared(tmp_path: Path) -> Path:
    """Create reproducible barrier inputs from the matched scan."""
    scan = tmp_path / "scan"
    _fixed_generator().generate(scan)
    preparer = _load(SOURCE / "prepare_barrier_benchmark.py", "barrier_preparer")
    output = tmp_path / "barriers"
    preparer.prepare(scan, output)
    return output


def test_prepare_barrier_inputs(prepared: Path) -> None:
    """Require nine mapped states, TS search, frequency, SP, and IRC inputs."""
    manifest = json.loads((prepared / "benchmark_manifest.json").read_text())
    assert manifest["purpose"].startswith("Unconstrained intrinsic gas-phase")
    assert manifest["states"] == ["reactant", "transition_state", "product"]
    assert len(manifest["records"]) == 9
    for record in manifest["records"]:
        calculation = prepared / record["directory"]
        gaussian = (calculation / "opt_freq_sp.gjf").read_text()
        assert "wB97XD/def2SVP" in gaussian
        assert "wB97XD/aug-cc-pVTZ SP" in gaussian
        assert "Freq" in gaussian
        assert len(set(record["reactive_atom_indices_1based"].values())) == 3
        if record["state"] == "transition_state":
            assert "Opt=(TS,CalcFC,NoEigenTest" in gaussian
            assert (calculation / "irc_forward.gjf").is_file()
            assert (calculation / "irc_reverse.gjf").is_file()
        else:
            assert "Opt=(CalcFC" in gaussian


def test_prepare_approximate_inputs(prepared: Path) -> None:
    """Convert FCHK geometries to matched XYZ and DFTB3 inputs."""
    converter = _load(
        SOURCE / "prepare_approximate_inputs.py", "prepare_approximate_inputs"
    )
    manifest = json.loads((prepared / "benchmark_manifest.json").read_text())
    for record in manifest["records"]:
        calculation = prepared / record["directory"]
        lines = (calculation / "initial.xyz").read_text().splitlines()
        atoms = [line.split() for line in lines[2:]]
        numbers = [{"H": 1, "C": 6, "N": 7, "O": 8}[atom[0]] for atom in atoms]
        coordinates = np.asarray(
            [[float(value) for value in atom[1:4]] for atom in atoms]
        )
        _write_fchk(calculation / "state.fchk", numbers, coordinates)
    converter.prepare(prepared)
    first = prepared / manifest["records"][0]["directory"]
    assert "DFT-optimized" in (first / "optimized.xyz").read_text()
    sqm = (first / "dftb3.in").read_text()
    assert "qm_theory='DFTB3'" in sqm
    assert f"qmcharge={manifest['records'][0]['charge']}" in sqm


def test_analyze_barriers_and_reject_bad_frequency(
    prepared: Path, tmp_path: Path
) -> None:
    """Validate IRC connectivity and matched forward-barrier arithmetic."""
    _load(SOURCE / "prepare_approximate_inputs.py", "prepare_approximate_inputs")
    analyzer = _load(SOURCE / "analyze_barrier_benchmark.py", "barrier_analyzer")
    manifest = json.loads((prepared / "benchmark_manifest.json").read_text())
    state_offsets = {"reactant": 0.0, "transition_state": 10.0, "product": 3.0}
    for record in manifest["records"]:
        calculation = prepared / record["directory"]
        offset = state_offsets[record["state"]]
        frequencies = (
            "-500.0 200.0 300.0"
            if record["state"] == "transition_state"
            else "100.0 200.0 300.0"
        )
        calculation.joinpath("opt_freq_sp.log").write_text(
            "Stationary point found.\n"
            f"Frequencies -- {frequencies}\n"
            "Zero-point correction= 0.100000\n"
            f"SCF Done: E(RwB97XD) = {-100.0 + offset / analyzer.HARTREE_TO_KCAL_MOL:.12f} A.U.\n"
            "Normal termination of Gaussian 16\n"
            f"SCF Done: E(RwB97XD) = {-200.0 + offset / analyzer.HARTREE_TO_KCAL_MOL:.12f} A.U.\n"
            "Normal termination of Gaussian 16\n"
        )
        calculation.joinpath("dftb3.out").write_text(
            f"Total SCF energy = {-1000.0 + 2.0 * offset:.8f} kcal/mol\n"
            "Calculation Completed\n"
        )
        calculation.joinpath("gfn2.out").write_text(
            f"| TOTAL ENERGY {-20.0 + 1.5 * offset / analyzer.HARTREE_TO_KCAL_MOL:.12f} Eh |\n"
            "normal termination of xtb\n"
        )
        lines = (calculation / "initial.xyz").read_text().splitlines()
        atoms = [line.split() for line in lines[2:]]
        numbers = [{"H": 1, "C": 6, "N": 7, "O": 8}[atom[0]] for atom in atoms]
        base = np.asarray([[float(value) for value in atom[1:4]] for atom in atoms])
        _write_fchk(calculation / "state.fchk", numbers, base)
        if record["state"] != "transition_state":
            continue
        indices = record["reactive_atom_indices_1based"]
        donor = base[indices["donor"] - 1]
        acceptor = base[indices["acceptor"] - 1]
        for direction, fraction in (("forward", 0.95 / 2.7), ("reverse", 1.75 / 2.7)):
            coordinates = base.copy()
            coordinates[indices["proton"] - 1] = donor + fraction * (acceptor - donor)
            _write_fchk(calculation / f"irc_{direction}.fchk", numbers, coordinates)
            calculation.joinpath(f"irc_{direction}.log").write_text(
                "Normal termination of Gaussian 16\n"
            )
    output = tmp_path / "analysis"
    analyzer.analyze(prepared, output)
    metrics = json.loads((output / "barrier_metrics.json").read_text())
    first = metrics["thr_to_terminal_amine"]
    assert np.isclose(
        first["methods"]["wb97xd_augccpvtz"]["forward_barrier_kcal_mol"], 10.0
    )
    assert np.isclose(first["barrier_errors_kcal_mol"]["dftb3"]["forward"], 10.0)

    bad = prepared / "thr_to_terminal_amine/reactant/opt_freq_sp.log"
    bad.write_text(bad.read_text().replace("100.0 200.0", "-100.0 200.0"))
    with pytest.raises(RuntimeError, match="significant imaginary"):
        analyzer.analyze(prepared, tmp_path / "rejected")


def test_reactive_mapping_rejects_nonmatching_coordinates() -> None:
    """Fail explicitly when a seed no longer matches the declared construction."""
    preparer = _load(SOURCE / "prepare_barrier_benchmark.py", "mapping_preparer")
    with pytest.raises(ValueError, match="Unable to recover"):
        preparer._reactive_indices(np.ones((3, 3)), 1.2)
