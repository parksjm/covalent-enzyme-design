"""Tests for preparing a constrained, crystal-referenced TI1 model."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from nyl12_covalent.ti1_bonded_model import build_ti1_bonded_model
from nyl12_covalent.ti1_reference_coordinates import (
    TORSION_LABELS,
    embed_bonded_reference,
)
from scripts import prepare_nyl12_ti1_bonded_reference as prepare


TARGETS = {"phi": -120.0, "psi": 130.0, "chi1": -60.0}


def _angular_error(actual: float, target: float) -> float:
    return abs((actual - target + 180.0) % 360.0 - 180.0)


def _embedded(cip: str = "R", seed: int = 721) -> tuple[Chem.Mol, dict[str, int]]:
    molecule, metadata = build_ti1_bonded_model(cip)
    atom_indices = metadata["atom_index_by_origin"]
    result = embed_bonded_reference(molecule, atom_indices, TARGETS, seed)
    return result, atom_indices


def test_bonded_embedding_preserves_cip_sets_torsions_and_is_deterministic() -> None:
    first, atom_indices = _embedded()
    second, _ = _embedded()

    assert first.GetNumAtoms() == 37
    assert first.GetNumConformers() == 1
    assert second.GetNumConformers() == 1
    for label, names in TORSION_LABELS.items():
        indices = [atom_indices[name] for name in names]
        angle = rdMolTransforms.GetDihedralDeg(first.GetConformer(), *indices)
        assert _angular_error(angle, TARGETS[label]) < 1e-5
    chiral_atom = first.GetAtomWithIdx(atom_indices["TIL:C*"])
    assert chiral_atom.GetProp("_CIPCode") == "R"
    np.testing.assert_allclose(
        np.asarray(first.GetConformer().GetPositions()),
        np.asarray(second.GetConformer().GetPositions()),
        rtol=0.0,
        atol=1e-8,
    )


def _reference(phi: float | None = -120.0) -> dict[str, object]:
    return {
        "source_path": "synthetic.cif",
        "sha256": "a" * 64,
        "chain": "A",
        "residue_number": 233,
        "phi_degrees": phi,
        "psi_degrees": 130.0,
        "chi1_degrees": -60.0,
        "previous_C_to_N_A": 1.34 if phi is not None else None,
        "C_to_next_N_A": 1.33,
        "thr_coordinates_A": {},
    }


def _prepare(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    output: Path,
    *,
    phi: float | None = -120.0,
    cores: int = 2,
    memory_gb: float = 4.0,
) -> None:
    cif = tmp_path / "source.cif"
    cif.write_text("synthetic CIF", encoding="utf-8")
    monkeypatch.setattr(
        prepare, "crystal_backbone_reference", lambda *_args: _reference(phi)
    )
    prepare.prepare_reference(cif, "A", 233, "R", None, 721, cores, memory_gb, output)


def test_prepare_writes_new_reference_and_gaussian_constraints(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    output = tmp_path / "reference-output"
    _prepare(monkeypatch, tmp_path, output)

    assert {path.name for path in output.iterdir()} == {
        "molecule.sdf",
        "coordinates.xyz",
        "atom_map.json",
        "input.gjf",
        "frequency.gjf",
        "esp.gjf",
    }
    metadata = json.loads((output / "atom_map.json").read_text(encoding="utf-8"))
    assert metadata["coordinate_reference"]["sha256"] == "a" * 64
    assert metadata["geometry_method"] == "wB97XD/def2TZVP"
    assert metadata["solvent"] == "gas phase"
    assert metadata["forcefield_validated"] is False
    assert metadata["full_stationary_minimum_claimed"] is False
    gjf = (output / "input.gjf").read_text(encoding="utf-8")
    assert "%nprocshared=2" in gjf.lower()
    assert "%mem=4" in gjf.lower()
    assert "wB97XD/def2TZVP" in gjf
    assert "Opt=(Tight,CalcFC,ModRedundant,MaxCycles=200)" in gjf
    assert "Freq" not in gjf
    assert "Freq" in (output / "frequency.gjf").read_text()
    assert "%chk=hessian.chk" in (output / "frequency.gjf").read_text()
    assert "HF/6-31G(d)" in (output / "esp.gjf").read_text()
    assert "%chk=esp.chk" in (output / "esp.gjf").read_text()
    constraints = [line for line in gjf.splitlines() if line.endswith(" F")]
    assert len(constraints) == 3
    atom_map = metadata["atom_index_by_origin"]
    for constraint, label in zip(constraints, TORSION_LABELS, strict=True):
        fields = constraint.split()
        expected_indices = [atom_map[name] + 1 for name in TORSION_LABELS[label]]
        assert [int(index) for index in fields[1:5]] == expected_indices
        assert fields[5] == "F"
        recorded = next(row for row in metadata["constraints"] if row["name"] == label)
        assert _angular_error(recorded["degrees"], TARGETS[label]) < 1e-5


def test_phi_is_required_or_explicit_cap_override_is_used(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    output = tmp_path / "missing-phi"
    cif = tmp_path / "source.cif"
    cif.write_text("synthetic CIF", encoding="utf-8")
    monkeypatch.setattr(
        prepare, "crystal_backbone_reference", lambda *_args: _reference(None)
    )
    with pytest.raises(ValueError):
        prepare.prepare_reference(cif, "A", 233, "R", None, 721, 2, 4.0, output)
    assert not output.exists()
    prepare.prepare_reference(cif, "A", 233, "R", -60.0, 721, 2, 4, output)
    metadata = json.loads((output / "atom_map.json").read_text())
    assert metadata["cap_phi_override_degrees"] == -60.0
    assert metadata["coordinate_reference"]["phi_degrees"] is None


@pytest.mark.parametrize(("cores", "memory"), [(0, 4.0), (2, 0.5)])
def test_invalid_resources_fail_before_creating_output(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    cores: int,
    memory: float,
) -> None:
    output = tmp_path / "invalid-resources"
    with pytest.raises(ValueError):
        _prepare(monkeypatch, tmp_path, output, cores=cores, memory_gb=memory)
    assert not output.exists()


def test_existing_output_is_never_overwritten(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    output = tmp_path / "existing"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("preserve", encoding="utf-8")
    with pytest.raises((FileExistsError, ValueError)):
        _prepare(monkeypatch, tmp_path, output)
    assert sentinel.read_text(encoding="utf-8") == "preserve"


def test_cli_applies_cap_override_and_writes_analysis_inputs(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    cif = tmp_path / "source.cif"
    cif.write_text("synthetic CIF", encoding="utf-8")
    output = tmp_path / "cli-reference"
    monkeypatch.setattr(
        prepare, "crystal_backbone_reference", lambda *_args: _reference(None)
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prepare_nyl12_ti1_bonded_reference",
            "--cif",
            str(cif),
            "--chain",
            "A",
            "--residue",
            "233",
            "--cip",
            "R",
            "--cap-phi",
            "-110",
            "--seed",
            "721",
            "--cores",
            "3",
            "--memory-gb",
            "5",
            "--output",
            str(output),
        ],
    )

    prepare.main()

    metadata = json.loads((output / "atom_map.json").read_text(encoding="utf-8"))
    assert metadata["cap_phi_override_degrees"] == -110.0
    assert metadata["coordinate_reference"]["sha256"] == "a" * 64
    assert metadata["cores"] == 3
    assert metadata["memory_gb"] == 5
    optimization = (output / "input.gjf").read_text(encoding="utf-8")
    assert "ModRedundant" in optimization
    assert sum(line.endswith(" F") for line in optimization.splitlines()) == 3

    esp = (output / "esp.gjf").read_text(encoding="utf-8")
    frequency = (output / "frequency.gjf").read_text(encoding="utf-8")
    assert "HF/6-31G(d)" in esp and "Pop=MK" in esp
    assert "wB97XD/def2TZVP" in frequency and "Freq" in frequency
