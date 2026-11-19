"""Tests for the small-model proton-transfer benchmark generator."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType

import numpy as np


class _FakeAxis:
    """Minimal Matplotlib axis used to exercise plotting logic."""

    def plot(self, *args: object, **kwargs: object) -> None:
        """Accept a plotted profile."""

    def axhline(self, *args: object, **kwargs: object) -> None:
        """Accept a reference line."""

    def set_title(self, title: str) -> None:
        """Accept an axis title."""

    def set_xlabel(self, label: str) -> None:
        """Accept an x-axis label."""

    def set_ylabel(self, label: str) -> None:
        """Accept a y-axis label."""

    def get_legend_handles_labels(self) -> tuple[list[object], list[str]]:
        """Return placeholder legend content."""
        return [], []


class _FakeFigure:
    """Minimal figure that records the requested output path."""

    def legend(self, *args: object, **kwargs: object) -> None:
        """Accept a figure legend."""

    def suptitle(self, *args: object, **kwargs: object) -> None:
        """Accept a figure title."""

    def tight_layout(self, *args: object, **kwargs: object) -> None:
        """Accept layout adjustment."""

    def savefig(self, path: Path, **kwargs: object) -> None:
        """Write a placeholder plot artifact."""
        path.write_bytes(b"plot")


def _install_fake_matplotlib() -> None:
    """Install a minimal pyplot module for dependency-free plot testing."""
    matplotlib = ModuleType("matplotlib")
    pyplot = ModuleType("matplotlib.pyplot")
    pyplot.subplots = lambda rows, columns, **kwargs: (
        _FakeFigure(),
        [[_FakeAxis() for _ in range(columns)] for _ in range(rows)],
    )
    pyplot.close = lambda figure: None
    matplotlib.pyplot = pyplot
    sys.modules["matplotlib"] = matplotlib
    sys.modules["matplotlib.pyplot"] = pyplot


def _load_script(name: str) -> ModuleType:
    """Load a scientific benchmark script without making it a package API."""
    path = (
        Path(__file__).parents[1]
        / f"qmmm/benchmarks/01_small_model_proton_transfer/source/{name}.py"
    )
    spec = importlib.util.spec_from_file_location(f"proton_benchmark_{name}", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load benchmark generator: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_generate_matched_benchmark_inputs(tmp_path: Path) -> None:
    """Generate all models and require complete matched method inputs."""
    module = _load_script("generate_benchmark")
    output = tmp_path / "generated"
    module.generate(output)

    manifest = json.loads((output / "benchmark_manifest.json").read_text())
    assert manifest["purpose"].startswith("Fixed-geometry")
    assert len(manifest["systems"]) == 3
    assert manifest["distances_A"] == [2.5, 2.7, 3.0]
    assert len(manifest["records"]) == 81

    for record in manifest["records"]:
        calculation = output / record["directory"]
        assert (calculation / "geometry.xyz").is_file()
        assert (calculation / "dftb3.in").is_file()
        assert (calculation / "wb97xd_augccpvtz.gjf").is_file()
        assert int((calculation / "charge.txt").read_text()) == record["charge"]
        gaussian = (calculation / "wb97xd_augccpvtz.gjf").read_text()
        assert "wB97XD/aug-cc-pVTZ" in gaussian
        assert f"{record['charge']} 1" in gaussian

    first = manifest["records"][0]
    assert first["donor_H_A"] == 0.95
    assert first["acceptor_H_A"] == 1.55
    last = manifest["records"][-1]
    assert last["donor_H_A"] == 1.98
    assert last["acceptor_H_A"] == 1.02


def test_rotation_handles_parallel_opposite_and_general_vectors() -> None:
    """Require the geometry helper to cover all rotation branches."""
    module = _load_script("generate_benchmark")
    source = np.array([1.0, 0.0, 0.0])
    for target in (
        np.array([1.0, 0.0, 0.0]),
        np.array([-1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
    ):
        rotation = module._rotation_between(source, target)
        assert np.allclose(rotation @ source, target, atol=1.0e-10)
        assert np.allclose(rotation.T @ rotation, np.eye(3), atol=1.0e-10)


def test_analyze_matched_outputs(tmp_path: Path) -> None:
    """Parse synthetic successful engine outputs and write comparison artifacts."""
    generator = _load_script("generate_benchmark")
    _install_fake_matplotlib()
    analyzer = _load_script("analyze_benchmark")
    generated = tmp_path / "generated"
    generator.generate(generated)
    manifest = json.loads((generated / "benchmark_manifest.json").read_text())
    for record in manifest["records"]:
        calculation = generated / record["directory"]
        point = record["point"]
        (calculation / "dftb3.out").write_text(
            f"Total SCF energy = {-1000.0 + 2.0 * point:.8f} kcal/mol\n"
            "Calculation Completed\n"
        )
        (calculation / "gfn2.out").write_text(
            f"| TOTAL ENERGY {-20.0 + 3.0 * point / analyzer.HARTREE_TO_KCAL_MOL:.12f} Eh |\n"
            "normal termination of xtb\n"
        )
        (calculation / "wb97xd_augccpvtz.log").write_text(
            f"SCF Done: E(RwB97XD) = {-30.0 + 2.5 * point / analyzer.HARTREE_TO_KCAL_MOL:.12f} A.U.\n"
            "Normal termination of Gaussian 16\n"
        )

    output = tmp_path / "analysis"
    analyzer.analyze(generated, output)
    metrics = json.loads((output / "profile_error_metrics.json").read_text())
    assert len(metrics) == 9
    assert (output / "matched_energies.csv").is_file()
    assert (output / "proton_transfer_method_comparison.png").is_file()
    first = metrics["thr_to_terminal_amine/R2.5"]
    assert np.isclose(
        first["dftb3"]["endpoint_reaction_energy_error_kcal_mol"], -4.0
    )
    assert np.isclose(
        first["gfn2_xtb"]["endpoint_reaction_energy_error_kcal_mol"], 4.0
    )
