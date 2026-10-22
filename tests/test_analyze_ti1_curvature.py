"""Tests for the reference TI1 constrained-curvature assessment."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from rdkit import Chem

from nyl12_covalent.ti1_bonded_model import build_ti1_bonded_model
from nyl12_covalent.ti1_reference_coordinates import (
    TORSION_LABELS,
    embed_bonded_reference,
)
from scripts import analyze_nyl12_ti1_curvature as analysis


TARGETS: dict[str, float] = {
    "phi": -120.0,
    "psi": 134.110276428,
    "chi1": -58.964417878,
}
INPUT_NAMES: tuple[str, ...] = (
    "reference.chk",
    "atom_map.json",
    "molecule.sdf",
    "input.gjf",
)
ATOMIC_UNITS_TO_ANGSTROM = 0.529177210903


def _sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _prepare_frequency_dir(directory: Path) -> tuple[Path, np.ndarray, np.ndarray]:
    """Write a self-consistent synthetic reference calculation input bundle."""
    directory.mkdir(parents=True)
    molecule, metadata = build_ti1_bonded_model("R")
    molecule = embed_bonded_reference(
        molecule,
        metadata["atom_index_by_origin"],
        TARGETS,
        20261022,
    )
    atom_map = dict(metadata)
    atom_map["geometry_method"] = "wB97XD/def2TZVP"
    atom_map["constraints"] = [
        {
            "name": name,
            "atoms_one_based": [
                atom_map["atom_index_by_origin"][label] + 1 for label in labels
            ],
            "degrees": TARGETS[name],
        }
        for name, labels in TORSION_LABELS.items()
    ]
    (directory / "atom_map.json").write_text(json.dumps(atom_map))
    writer = Chem.SDWriter(str(directory / "molecule.sdf"))
    writer.write(molecule)
    writer.close()
    (directory / "input.gjf").write_text("Synthetic Gaussian input for test.\n")
    (directory / "reference.chk").write_bytes(b"synthetic checkpoint placeholder")
    (directory / "hessian.fchk").write_text(
        "synthetic formatted checkpoint placeholder\n"
    )
    (directory / "scientific_acceptance.json").write_text(
        json.dumps(
            {
                "accepted": True,
                "checkpoint_geometry_match": {"accepted": True},
            }
        )
    )
    (directory / "process_result.json").write_text(json.dumps({"exit_status": 0}))
    identity = {name: _sha256(directory / name) for name in INPUT_NAMES}
    (directory / "input_identity.json").write_text(json.dumps(identity))

    conformer = molecule.GetConformer()
    coordinates_bohr = np.asarray(conformer.GetPositions()) / ATOMIC_UNITS_TO_ANGSTROM
    atomic_numbers = np.asarray(
        [atom.GetAtomicNum() for atom in molecule.GetAtoms()], dtype=int
    )
    return directory, atomic_numbers, coordinates_bohr


def _install_synthetic_hessian(
    monkeypatch: pytest.MonkeyPatch,
    atomic_numbers: np.ndarray,
    coordinates_bohr: np.ndarray,
    *,
    gradient: np.ndarray | None = None,
) -> None:
    """Supply the Hessian reader with a known positive quadratic model."""
    atom_count = len(atomic_numbers)
    if gradient is None:
        gradient = np.zeros((atom_count, 3))
    monkeypatch.setattr(
        analysis,
        "read_gaussian_fchk",
        lambda path, **identity_checks: SimpleNamespace(
            atomic_numbers=atomic_numbers,
            coordinates_bohr=coordinates_bohr,
            hessian_eh_per_bohr2=np.eye(3 * atom_count),
            gradient_eh_per_bohr=gradient,
        ),
    )


def _run_analysis(
    frequency_dir: Path,
    output: Path,
    *,
    step_bohr: float = 2e-4,
    stationarity_tolerance: float = 3e-5,
    curvature_tolerance: float = 1e-6,
    derivative_relative_tolerance: float = 1e-3,
) -> dict[str, Any]:
    """Run the analyzer with explicit numerical settings."""
    return analysis.analyze_reference_curvature(
        frequency_dir,
        output,
        step_bohr=step_bohr,
        stationarity_tolerance=stationarity_tolerance,
        curvature_tolerance=curvature_tolerance,
        derivative_relative_tolerance=derivative_relative_tolerance,
    )


def test_synthetic_quadratic_reference_is_accepted_in_stated_scope(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A manufactured quadratic model validates plumbing, not real QM data."""
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(monkeypatch, atomic_numbers, coordinates)

    result = _run_analysis(frequency_dir, tmp_path / "assessment")

    assert result["accepted"] is True
    assert result["scope"].startswith("Constrained reference stability only")
    assert result["forcefield_validated"] is False
    assert result["fine"]["tangent_dimension"] == 102
    assert (tmp_path / "assessment" / "assessment.json").is_file()


def test_missing_gradient_fails_before_creating_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The analyzer requires a gradient for its stationarity assessment."""
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(
        monkeypatch,
        atomic_numbers,
        coordinates,
        gradient=np.full((len(atomic_numbers), 3), np.nan),
    )
    # A missing gradient is represented by the parser's explicit None value.
    monkeypatch.setattr(
        analysis,
        "read_gaussian_fchk",
        lambda path, **identity_checks: SimpleNamespace(
            atomic_numbers=atomic_numbers,
            coordinates_bohr=coordinates,
            hessian_eh_per_bohr2=np.eye(3 * len(atomic_numbers)),
            gradient_eh_per_bohr=None,
        ),
    )
    output = tmp_path / "assessment"

    with pytest.raises(ValueError):
        _run_analysis(frequency_dir, output)

    assert not output.exists()


@pytest.mark.parametrize(
    ("manifest", "filename", "replacement"),
    [
        (
            "scientific_acceptance.json",
            None,
            {"accepted": False, "checkpoint_geometry_match": {"accepted": True}},
        ),
        ("process_result.json", None, {"exit_status": 1}),
    ],
)
def test_unaccepted_or_failed_reference_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    manifest: str,
    filename: str | None,
    replacement: dict[str, Any],
) -> None:
    """Raw scientific acceptance and process success gate analysis."""
    del filename
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(monkeypatch, atomic_numbers, coordinates)
    (frequency_dir / manifest).write_text(json.dumps(replacement))

    with pytest.raises(ValueError):
        _run_analysis(frequency_dir, tmp_path / "assessment")


def test_missing_required_acceptance_key_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Acceptance records must explicitly identify geometry agreement."""
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(monkeypatch, atomic_numbers, coordinates)
    (frequency_dir / "scientific_acceptance.json").write_text(
        json.dumps({"accepted": True})
    )

    with pytest.raises((KeyError, ValueError)):
        _run_analysis(frequency_dir, tmp_path / "assessment")


def test_changed_identity_input_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every hashed source input must still match its recorded identity."""
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(monkeypatch, atomic_numbers, coordinates)
    (frequency_dir / "input.gjf").write_text("changed after identity creation\n")

    with pytest.raises(ValueError):
        _run_analysis(frequency_dir, tmp_path / "assessment")


def test_invalid_constraint_atoms_are_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A malformed torsion definition cannot pass identity validation."""
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(monkeypatch, atomic_numbers, coordinates)
    metadata_path = frequency_dir / "atom_map.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["constraints"][0]["atoms_one_based"] = [1, 1, 2, 3]
    metadata_path.write_text(json.dumps(metadata))
    identity = json.loads((frequency_dir / "input_identity.json").read_text())
    identity["atom_map.json"] = _sha256(metadata_path)
    (frequency_dir / "input_identity.json").write_text(json.dumps(identity))

    with pytest.raises(ValueError):
        _run_analysis(frequency_dir, tmp_path / "assessment")


def test_existing_output_directory_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The analyzer refuses to overwrite an existing assessment directory."""
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(monkeypatch, atomic_numbers, coordinates)
    output = tmp_path / "assessment"
    output.mkdir()

    with pytest.raises((FileExistsError, ValueError)):
        _run_analysis(frequency_dir, output)


def test_nonpositive_derivative_tolerance_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Finite-difference derivative acceptance requires positive tolerance."""
    frequency_dir, atomic_numbers, coordinates = _prepare_frequency_dir(
        tmp_path / "frequency"
    )
    _install_synthetic_hessian(monkeypatch, atomic_numbers, coordinates)

    with pytest.raises(ValueError):
        _run_analysis(
            frequency_dir,
            tmp_path / "assessment",
            derivative_relative_tolerance=0.0,
        )


def test_cli_returns_nonzero_for_rejected_assessment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """CLI propagates a rejected scientific assessment as exit status 3."""
    called: dict[str, Any] = {}

    def fake_analyze(
        frequency_dir: Path,
        output: Path,
        *,
        step_bohr: float,
        stationarity_tolerance: float,
        curvature_tolerance: float,
        derivative_relative_tolerance: float,
    ) -> dict[str, Any]:
        called.update(
            frequency_dir=frequency_dir,
            output=output,
            step_bohr=step_bohr,
            stationarity_tolerance=stationarity_tolerance,
            curvature_tolerance=curvature_tolerance,
            derivative_relative_tolerance=derivative_relative_tolerance,
        )
        return {"accepted": False}

    monkeypatch.setattr(analysis, "analyze_reference_curvature", fake_analyze)
    monkeypatch.setattr(
        "sys.argv",
        [
            "analyze_nyl12_ti1_curvature",
            "--frequency-dir",
            "freq",
            "--output",
            "out",
            "--step-bohr",
            "0.0002",
            "--stationarity-tolerance",
            "0.00003",
            "--curvature-tolerance",
            "0.000001",
            "--derivative-relative-tolerance",
            "0.001",
        ],
    )

    with pytest.raises(SystemExit) as error:
        analysis.main()

    assert error.value.code == 3
    assert called == {
        "frequency_dir": Path("freq"),
        "output": Path("out"),
        "step_bohr": 2e-4,
        "stationarity_tolerance": 3e-5,
        "curvature_tolerance": 1e-6,
        "derivative_relative_tolerance": 1e-3,
    }
