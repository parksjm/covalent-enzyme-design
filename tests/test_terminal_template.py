"""Tests for the neutral threonine methylamide template."""

import json

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from nyl12_covalent.terminal_template import (
    _audit_coordinate_stereochemistry,
    _audit_graph_stereochemistry,
    _build_template,
    _gaussian_input,
    _name_hydrogens,
    _write_atom_metadata,
    _write_xyz,
    prepare_templates,
)


def _embedded_template() -> tuple[Chem.Mol, dict[str, int]]:
    molecule, atom_indices = _build_template()
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 314159
    parameters.numThreads = 1
    assert AllChem.EmbedMolecule(molecule, parameters) == 0
    return molecule, atom_indices


def test_template_graph_has_expected_formula_charge_and_named_protons() -> None:
    molecule, atom_indices = _build_template()

    assert CalcMolFormula(molecule) == "C5H12N2O2"
    assert Chem.GetFormalCharge(molecule) == 0
    _audit_graph_stereochemistry(molecule, atom_indices)
    names = [atom.GetProp("atom_name") for atom in molecule.GetAtoms()]
    residues = [atom.GetProp("residue_name") for atom in molecule.GetAtoms()]
    assert names.count("H1") == 1
    assert names.count("H2") == 1
    assert names.count("HG21") == 1
    assert names.count("HH31") == 1
    nme_n = next(
        atom for atom in molecule.GetAtoms() if atom_indices["NME_N"] == atom.GetIdx()
    )
    assert nme_n.GetProp("atom_name") == "N"
    assert nme_n.GetProp("residue_name") == "NME"
    assert len(residues) == molecule.GetNumAtoms()


def test_hydrogen_naming_rejects_incomplete_inventory() -> None:
    molecule, atom_indices = _build_template()
    editable = Chem.RWMol(molecule)
    hydrogen_index = next(
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetSymbol() == "H" and atom.GetProp("atom_name") == "H1"
    )
    editable.RemoveAtom(hydrogen_index)
    with pytest.raises(ValueError, match="Hydrogen inventory mismatch"):
        _name_hydrogens(editable.GetMol(), atom_indices)


def test_coordinate_stereo_and_gaussian_routes_are_explicit() -> None:
    molecule, atom_indices = _embedded_template()
    conformer = molecule.GetConformer()
    rdMolTransforms.SetDihedralDeg(
        conformer,
        atom_indices["N"],
        atom_indices["CA"],
        atom_indices["C"],
        atom_indices["NME_N"],
        -40.0,
    )
    _audit_coordinate_stereochemistry(molecule, atom_indices)
    gaussian_input = _gaussian_input(molecule, "neutral_thr_nme.chk", "wb97xd-def2tzvp")

    assert "wB97XD/def2TZVP Int=UltraFine" in gaussian_input
    assert "EmpiricalDispersion" not in gaussian_input
    assert "--Link1--" in gaussian_input
    assert "HF/6-31G(d) Geom=AllCheck Guess=Read" in gaussian_input
    assert "Pop=MK IOp(6/33=2,6/42=6)" in gaussian_input
    assert "0 1" in gaussian_input
    b3lyp_input = _gaussian_input(
        molecule, "neutral_thr_nme.chk", "b3lyp-d3bj-def2tzvp"
    )
    assert "B3LYP/def2TZVP EmpiricalDispersion=GD3BJ" in b3lyp_input


def test_metadata_and_xyz_preserve_atom_order(tmp_path) -> None:
    molecule, _ = _embedded_template()
    atom_map_path = tmp_path / "atom_map.json"
    xyz_path = tmp_path / "coordinates.xyz"
    _write_atom_metadata(molecule, atom_map_path)
    _write_xyz(molecule, xyz_path)

    atom_map = json.loads(atom_map_path.read_text(encoding="utf-8"))
    xyz_lines = xyz_path.read_text(encoding="utf-8").splitlines()
    assert len(atom_map["atoms"]) == molecule.GetNumAtoms()
    assert int(xyz_lines[0]) == molecule.GetNumAtoms()
    assert len(xyz_lines[2:]) == molecule.GetNumAtoms()


def test_prepare_rejects_existing_directory_before_embedding(tmp_path) -> None:
    output = tmp_path / "existing"
    output.mkdir()

    with pytest.raises(FileExistsError):
        prepare_templates(output, "b3lyp-d3bj-def2tzvp")


def test_prepare_requires_supported_geometry_method(tmp_path) -> None:
    with pytest.raises(ValueError, match="Unsupported geometry method"):
        prepare_templates(tmp_path / "new", "hf-6-31g-star")  # type: ignore[arg-type]


def test_prepare_writes_seeded_conformer_artifacts(tmp_path) -> None:
    output = prepare_templates(
        tmp_path / "prepared",
        geometry_method="wb97xd-def2tzvp",
        seeds=(271828,),
    )

    conformer_dir = output / "conformer_01"
    assert sorted(path.name for path in conformer_dir.iterdir()) == [
        "atom_map.json",
        "coordinates.xyz",
        "input.gjf",
        "molecule.sdf",
    ]
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["formula"] == "C5H12N2O2"
    assert manifest["formal_charge"] == 0
    assert manifest["conformers"][0]["seed"] == 271828
    assert manifest["conformers"][0]["psi_degrees"] == -40.0
